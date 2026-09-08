#仅供软著申请用！
#合并所有源代码文件
import os
target_file_patterns = [
    {"dir": "setup.py", "extensions": [".py"]},
    {"dir": "scripts", "extensions": [".sh"]},
    {"dir": "plugin_packs", "extensions": [".json", ".py", ".sh"]},
    {"dir": "ipamd", "extensions": [".py", ".xml"]},
    {"dir": "demo", "extensions": [".py"]},
]
output_file_name = 'software_copyright_material/src.txt'
output_file = open(output_file_name, 'w', encoding='utf-8')
comment_char = '###'

def process_file(file, extensions):
    line_cnt = 0
    is_dir = os.path.isdir(file)
    if is_dir:
        files = os.listdir(file)
        for sub_file in files:
            path = os.path.join(file, sub_file)
            line_cnt += process_file(path, extensions)
    else:
        _, extension = os.path.splitext(file)
        if extension in extensions:
            with open(file, 'r', encoding='utf-8') as f:
                content = f.read().strip()
            line_cnt += len(content.split('\n'))
            output_file.write(f'{comment_char} File: {file}\n')
            output_file.write(content + '\n\n')
    return line_cnt

total_lines = 0
for target_file_pattern in target_file_patterns:
    target_dir = target_file_pattern['dir']
    target_extensions = target_file_pattern['extensions']
    total_lines += process_file(target_dir, target_extensions)
print(f'Total {total_lines} lines.')
print(f'Src file {output_file_name} generated.')
