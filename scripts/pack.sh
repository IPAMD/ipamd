#!/bin/bash

# 列出所有已安装的插件
plugins=$(python ipamd/cli.py -l)

# 对每个插件执行卸载命令
echo "Uninstalling all plugins..."
while IFS= read -r plugin; do
    if [ -n "$plugin" ]; then
        echo "Uninstalling plugin: $plugin"
        python ipamd/cli.py -r "$plugin"
    fi
done <<< "$plugins"

echo "All plugins have been uninstalled."

# 进入plugin_packs目录并打包所有插件
cd plugin_packs

# 创建zipped目录（如果不存在）
mkdir -p zipped

# 对除了zipped的文件夹运行打包命令
echo "Packing all plugin directories..."
for dir in */; do
    # 移除目录末尾的斜杠
    dir_name=$(basename "$dir")
    
    # 跳过zipped目录
    if [ "$dir_name" != "zipped" ]; then
        echo "Packing plugin: $dir_name"
        python ../ipamd/cli.py -p "$dir_name"
    fi
done

# 将生成的zip文件移动到zipped文件夹
echo "Moving zip files to zipped directory..."
mv *.zip zipped/ 2>/dev/null

echo "All plugins have been packed and moved to zipped directory."

cd ..
python setup.py sdist