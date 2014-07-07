# -*- python -*-
import os, sys

from lsst.sconsUtils import dependencies 
from lsst.sconsUtils import env

# Command line options
vars = Variables("custom.py")
vars.AddVariables(PathVariable('prefix', 'Installation prefix', '.', PathVariable.PathIsDirCreate),
                  BoolVariable('debug', 'Debug build?', 0)
                 )

dependencies.configure("partition")
print "XXXXXXXXX %s" % env.getLibs("main")
print "XXXXXXXXX %s" % env.subst('$LIBS')
print "XXXXXXXXX LIBPATH : %s" % env.subst('$LIBPATH')
env.Append(LIBS=[env.getLibs("main")])

# Build environment setup
#env = Environment(variables=vars)
env.Help("Spherical data partitioning and duplication utilities.")

env.Append(CCFLAGS=['-g', '-Wall', '-Wextra'])
if env['debug']:
    variant_dir = 'build/debug'
else:
    variant_dir = 'build/release'
    env.Append(CCFLAGS=['-O2'])
if not os.path.isabs(env['prefix']):
    # `scons prefix=<dir> install` seems buggy when <dir> is relative
    env['prefix'] = os.path.normpath(os.path.join(GetLaunchDir(), env['prefix']))

env.Alias('install', ['$prefix/bin', '$prefix/ups'])

# Recurse
if not GetOption('help'):
    SConscript(dirs='.', variant_dir=variant_dir, exports='env')

