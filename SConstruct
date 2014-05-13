# -*- python -*-
import os, sys

from lsst.sconsUtils import dependencies
from lsst.sconsUtils import env

dependencies.configure("partition")

env.Help("Spherical data partitioning and duplication utilities.")

env.Append(CCFLAGS=['-Wextra'])
if env['debug']:
    # -g flag is handled by sconsUtils
    build_cfg = 'debug'
else:
    build_cfg = 'release'
    env.Append(CCFLAGS=['-O2'])

if not os.path.isabs(env['prefix']):
    # `scons prefix=<dir> install` seems buggy when <dir> is relative
    env['prefix'] = os.path.normpath(os.path.join(GetLaunchDir(), env['prefix']))

env.Alias('install', ['$prefix/bin'])

# Recurse
if not GetOption('help'):
    for dir in ['src', 'tests']:
        variant_dir = os.path.join("build", build_cfg, dir)
        SConscript(dirs=[dir], variant_dir=variant_dir)

