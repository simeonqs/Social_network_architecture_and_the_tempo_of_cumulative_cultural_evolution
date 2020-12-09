This is the Python code to run the agent based models. 

To run, navigate to containing directory and run using python3, takes one argument which is the graph type to use
graph types: clustered, degree, full, modular_and_clustered, modular, multilevel, small_world

Example:
```bash
cd "/YOUR/PATH/code/2_Python_agent_based_models/code"
python3 M1.py "clustered"
python3 M2.py "smallworld"
```
This code requires the following modules:

* numpy 
* networkx 
* random
* warnings
* os
* sys

You can install them with:
```bash
pip3 install numpy
```

Latest tested on:

{'commit_hash': 'f7f2eae63',
 'commit_source': 'installation',
 'default_encoding': 'UTF-8',
 'ipython_path': '/usr/local/lib/python3.7/site-packages/IPython',
 'ipython_version': '7.10.2',
 'os_name': 'posix',
 'platform': 'Darwin-18.7.0-x86_64-i386-64bit',
 'sys_executable': '/usr/local/opt/python/bin/python3.7',
 'sys_platform': 'darwin',
 'sys_version': '3.7.7 (default, Mar 10 2020, 15:43:03) \n'
                '[Clang 11.0.0 (clang-1100.0.33.17)]'}
