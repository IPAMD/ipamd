import copy
import json
from multiprocessing import Process, Pipe
import re
import signal
import os
import time
from importlib import import_module
from rich.progress import Progress, BarColumn, TimeElapsedColumn, TimeRemainingColumn
from numba import cuda
from ipamd.public.utils.parser import value_of
from ipamd.public.utils.plugin_manager import PluginBase
from ipamd.public.utils.output import error, captured, info, captured_async, verbose
from ipamd.public import shared_data


class Simulation:
    class __Parser(PluginBase):
        def __init__(self):
            extra_plugin_dir = shared_data.config.get('simulation_plugin_dir')
            plugin_dir = ([
                shared_data.module_installation_dir + '/plugins/integrator',
                shared_data.module_installation_dir + '/plugins/force'
            ] + extra_plugin_dir)
            super().__init__(plugin_dir)
            self.load_all()

    def __init__(self, app):
        self.__app = app
        self.__parser = self.__Parser()
        @captured
        def f():
            return import_module('ipamd.app.lib.galamost')

        self.__gala_core, _ = f()

    @property
    def working_dir(self):
        """
        The working directory of the simulation.
        """
        return self.__app.working_dir

    def new_simulation(self,
                       simulation_box,
                       job_name='simulation',
                       dt=0.01,
                       total_time=0,
                       snap_shot=0,
                       run_step=0, period=0,
                       thermo_bath='langevin_nvt',
                       res_auto_read=True,
                       minimize_energy=True,
                       fixed_particle=[]
        ):
        """
        Create a new simulation job.

        @param simulation_box: The simulation box.
        @param job_name: The name of the simulation job.
        @param dt: The time step of the simulation.
        @param total_time: The total time of the simulation.
        @param snap_shot: The snapshot interval of the simulation.
        @param run_step: The total number of steps of the simulation.
        @param period: The period of the simulation.
        @param thermo_bath: The thermal bath of the simulation.
        @param res_auto_read: Whether to automatically read the results of the simulation.
        @param minimize_energy: Whether to minimize the energy of the simulation.
        """
        if total_time != 0 and snap_shot!=0:
            run_step = int(total_time / dt)
            period = int(run_step / snap_shot)

        return self.__Simulation_Job(
            simulation_box=simulation_box,
            job_name=job_name,
            run_step=run_step,
            dt=dt,
            period=period,
            parser=self.__parser,
            working_dir=self.working_dir,
            gala_core=self.__gala_core,
            gpu_id=self.__app.gpu_id,
            thermo_bath=thermo_bath,
            force_field=self.__app.force_field,
            minimize_energy=minimize_energy,
            fixed_particle=fixed_particle
        )

    class __Simulation_Job:
        """
        Simulation Job Class
        """
        def __init__(
            self,
            simulation_box,
            job_name,
            run_step,
            dt,
            parser,
            period,
            working_dir,
            gala_core,
            gpu_id,
            thermo_bath,
            force_field,
            minimize_energy,
            fixed_particle
        ):
            self.__total_step = run_step
            self.__target_step = run_step
            self.__gpu_id = 0
            self.__dt = dt
            self.__simulation_box = simulation_box
            self.job_name = job_name
            self.__parser = parser
            self.__working_dir = working_dir
            self.__gala_core = gala_core
            self.__gpu_id = gpu_id
            self.period = period
            self.__thermo_bath = thermo_bath
            self.__force_field = force_field
            self.__should_minimize_energy = minimize_energy
            self.__fixed_particle = fixed_particle

        def __parse_force(self, all_info, force, param):
            parser = getattr(self.__parser, force)
            return parser(param, all_info, self.__gala_core)

        def __parse_integrator(self, all_info, group, integrator, param):
            parser = getattr(self.__parser, integrator)
            return parser(param, all_info, group, self.__gala_core)

        def __outputs(self):
            files = os.listdir(self.__working_dir)
            pattern = '^' + self.job_name + r'\.\d+\.xml$'
            file_list = [f for f in files if re.match(pattern, f)]
            sorted_list = sorted(file_list, key=lambda x: int(x.split('.')[-2]))
            return sorted_list

        def __latest(self, file_list):
            max_step = -1
            latest_output = self.job_name
            for file in file_list:
                digit = int(re.findall(r'\d+', file)[-1])
                if digit > max_step:
                    max_step = digit
                    latest_output = file[:-4]
            return max_step, latest_output

        def clean_cache(self, remove_info=True):
            """
            Clean the data of the previous simulation.
            @param remove_info: Whether to remove the info file.
            """
            matching_files = self.__outputs()
            for file in matching_files:
                os.remove(os.path.join(self.__working_dir, file))
            log_filename = os.path.join(self.__working_dir, self.job_name + '.log')
            if os.path.exists(log_filename):
                os.remove(log_filename)
            if os.path.exists(os.path.join(self.__working_dir, self.job_name + '.xml')):
                os.remove(os.path.join(self.__working_dir, self.job_name + '.xml'))
            if os.path.exists(os.path.join(self.__working_dir, self.job_name + '.map')):
                os.remove(os.path.join(self.__working_dir, self.job_name + '.map'))
            if remove_info and os.path.exists(
                os.path.join(self.__working_dir, self.job_name + '.info')
            ):
                os.remove(os.path.join(self.__working_dir, self.job_name + '.info'))

        def exist(self):
            """
            Check if the simulation exists.
            @return: Whether the simulation exists.
            """
            info_file_path = os.path.join(self.__working_dir, self.job_name + '.info')
            return os.path.exists(info_file_path)

        def run(self, yes=''):
            """
            Run the simulation.
            @param yes: skip the interaction with the user.
            """
            ready_to_run = True
            input_file_path = os.path.join(self.__working_dir, self.job_name + '.xml')
            new_simulation = False
            info_file_path = os.path.join(self.__working_dir, self.job_name + '.info')
            if self.exist():
                with open(info_file_path, 'r', encoding='utf-8') as f:
                    simulation_info = json.load(f)
                    run_step = simulation_info['run_step']
                    matching_files = self.__outputs()
                    max_step, latest_output = self.__latest(matching_files)
                    def make_sure(hint, yes, no, override):
                        if override:
                            return override == yes
                        return input(f"{hint} ({yes}/{no}): ")  == yes

                    if max_step < run_step and max_step != -1:
                        if make_sure('Previous simulation not finished, continue?', 'y', 'n', yes):
                            input_file_path = os.path.join(self.__working_dir, latest_output + '.xml')
                            remaining_step = run_step - max_step
                            self.__total_step = remaining_step
                            self.__target_step = run_step
                        else:
                            self.clean_cache()
                            new_simulation = True
                    else:
                        if make_sure('Find previous simulation, rerun?', 'y', 'n', yes):
                            self.clean_cache()
                            new_simulation = True
                        else:
                            ready_to_run = False
            else:
                new_simulation = True

            if new_simulation:
                with open(info_file_path, 'w', encoding='utf-8') as f:
                    data = {
                        'job_name': self.job_name,
                        'run_step': self.__total_step,
                        'dt': self.__dt,
                        'period': self.period
                    }
                    json.dump(data, f, indent=4)

            if ready_to_run:
                cuda.close()
                self.__simulation_box.to_xml(self.job_name)
                if new_simulation and self.__should_minimize_energy:
                    self.__minimize_energy(input_file_path)
                    em_main_filename = f'{self.job_name}_em'
                    em_file_path = os.path.join(self.__working_dir, f'{em_main_filename}.xml')
                    info(f'writing {em_file_path} with energy minimized configuration')
                    matching_files = self.__outputs()
                    max_step, latest_output = self.__latest(matching_files)
                    self.__simulation_box.clean()
                    self.__simulation_box.new_frame()
                    self.__simulation_box.read_xml(latest_output)
                    self.__simulation_box.to_xml(em_main_filename)
                    input_file_path = em_file_path
                    for old_file in matching_files:
                        os.remove(os.path.join(self.__working_dir, old_file))
                self.__run_simulation(input_file_path)

            self.__simulation_box.clean()
            matching_files = self.__outputs()
            for file in matching_files:
                self.__simulation_box.new_frame()
                self.__simulation_box.read_xml(file[:-4])

            info('Simulation finished')

        @captured
        def __init_info(self, path):
            build_method = self.__gala_core.XMLReader(path)
            perform_config = self.__gala_core.PerformConfig(str(self.__gpu_id))
            return self.__gala_core.AllInfo(build_method, perform_config)

        @captured
        def __init_app(self, all_info, dt):
            return self.__gala_core.Application(all_info, dt)

        @captured
        def __init_force(self, all_info, target_app):
            lookup_table = self.__simulation_box.env.values
            for force in self.__force_field.ff_param:
                params = copy.deepcopy(self.__force_field.ff_param[force])
                for key in params.keys():
                    if isinstance(params[key], dict):
                        for sub_key in params[key].keys():
                            params[key][sub_key] = value_of(params[key][sub_key], lookup_table)
                    else:
                        params[key] = value_of(params[key], lookup_table)
                f = self.__parse_force(all_info, force, params)
                if isinstance(f, list):
                    for sf in f:
                        target_app.add(sf)
                else:
                    target_app.add(f)

            return 0 #in order to avoid warning

        @captured
        def __init_essentials(self, all_info, target_app):
            sort_method = self.__gala_core.Sort(all_info)
            sort_method.setPeriod(500)
            target_app.add(sort_method)

            zero_momentum = self.__gala_core.ZeroMomentum(all_info)
            zero_momentum.setPeriod(100)
            target_app.add(zero_momentum)
            return 0  # in order to avoid warning

        @captured
        def __init_log(self, all_info, target_app, period):
            group = self.__gala_core.ParticleSet(all_info, "all")
            comp_info = self.__gala_core.ComputeInfo(all_info, group)
            d_info = self.__gala_core.DumpInfo(
                all_info,
                comp_info,
                self.__working_dir + '/' + self.job_name + '.log'
            )
            d_info.setPeriod(int(period))
            target_app.add(d_info)

        @captured
        def __init_snapshot(self, all_info, target_app, period):
            xml = self.__gala_core.XMLDump(all_info, self.__working_dir + '/' + self.job_name)
            xml.setPeriod(int(period))
            xml.setOutput(
                [
                    'bond',
                    'image',
                    'molecule',
                    'velocity',
                    'mass',
                    'charge',
                    'body',
                    'angle'
                ]
            )
            target_app.add(xml)

        @captured
        def __init_integrator(self, all_info, target_app, type_='langevin_nvt'):
            rigid_body = False
            prop = self.__simulation_box.current_frame().properties()
            for molecule in prop['molecules']:
                for rigid in molecule['rigid_group']:
                    if rigid != -1:
                        rigid_body = True
                        break
            exclude_list = [-index for index in self.__fixed_particle]
            if rigid_body:
                group_rigid = self.__gala_core.ParticleSet(all_info, ['body'] + exclude_list)
                integrator_r = self.__parse_integrator(all_info, group_rigid, type_, {
                    'temperature': self.__simulation_box.env.values['temperature'],
                    'pressure': self.__simulation_box.env.values['pressure'],
                    'rigid': True
                })
                target_app.add(integrator_r)
            group_non_rigid = self.__gala_core.ParticleSet(all_info, ['non_body'] + exclude_list)
            integrator = self.__parse_integrator(all_info, group_non_rigid, type_, {
                'temperature': self.__simulation_box.env.values['temperature'],
                'pressure': self.__simulation_box.env.values['pressure'],
                'rigid': False
            })
            target_app.add(integrator)

        def __minimize_energy(self, input_file_path):
            verbose('reading input file ' + input_file_path)
            em_step = 1e4
            all_info, _ = self.__init_info(input_file_path)
            app, _ = self.__init_app(all_info, dt=0.001)
            self.__init_force(all_info, app)
            self.__init_integrator(all_info, app, type_='em')
            self.__init_essentials(all_info, app)
            self.__init_snapshot(all_info, app, em_step)
            verbose('minimizing energy')

            @captured
            def run_md():
                app.run(int(em_step))
            run_md()

        def __run_simulation(self, input_file_path):
            all_info, _ = self.__init_info(input_file_path)
            app, _ = self.__init_app(all_info, dt=self.__dt)
            self.__init_force(all_info, app)
            self.__init_essentials(all_info, app)
            self.__init_log(all_info, app, self.period)
            self.__init_snapshot(all_info, app, self.period)
            self.__init_integrator(all_info, app, type_=self.__thermo_bath)

            parent_conn, child_conn = Pipe()

            printer_process = Process(target=self.__printer, args=(child_conn,))
            printer_process.start()

            @captured_async(parent_conn, timeout=0.1)
            def run_md():
                app.run(self.__total_step)
            run_md()
            printer_process.join()

        def __printer(self, conn):
            def exit_(signum, frame):
                os.kill(os.getppid(), signal.SIGKILL)
                os._exit(0)

            signal.signal(signal.SIGINT, exit_)
            signal.signal(signal.SIGTERM, exit_)


            with Progress(
                    '[progress.description]{task.description}',
                    BarColumn(),
                    "[progress.percentage]{task.percentage:>3.2f}%",
                    "[yellow]⏱",
                    TimeElapsedColumn(),
                    "[cyan]⏳",
                    TimeRemainingColumn(),
                    refresh_per_second=1) as progress:

                task = None
                simulation_finished = False
                while not simulation_finished:
                    message = conn.recv()
                    for line in message.split('\n'):
                        if '--start--' in line:
                            task = progress.add_task(self.job_name, total=self.__target_step)
                            progress.update(task,
                                            description=f'{self.job_name}(TPS: *)',
                                            completed=200)
                        elif 'STATUS RUNNING' in line:
                            numbers = list(map(str, re.findall(r'\d+\.?\d*', line)))
                            tps = float(numbers[0])
                            current = int(numbers[1])
                            progress.update(task,
                                            description=f'{self.job_name}(TPS: {int(round(tps))})',
                                            completed=current)
                        elif '--end--' in line:
                            progress.update(task, completed=self.__total_step)
                            simulation_finished = True
                        elif 'Error' in line:
                            simulation_finished = True
                            error(f'{line[10:]}')
                    time.sleep(0.1)
                conn.close()
