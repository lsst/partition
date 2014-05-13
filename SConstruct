# -*- python -*-
import os, sys

# Command line options
BOOST_DIR = os.environ.get('BOOST_DIR')
BOOST_INC = None
BOOST_LIB = None
if BOOST_DIR:
    BOOST_INC = os.path.join(BOOST_DIR, 'include')
    BOOST_LIB = os.path.join(BOOST_DIR, 'lib')
vars = Variables("custom.py")
vars.AddVariables(PathVariable('prefix', 'Installation prefix', '.', PathVariable.PathIsDirCreate),
                  BoolVariable('debug', 'Debug build?', 0),
                  PathVariable('BOOST_INC', 'Boost include directory', BOOST_INC, PathVariable.PathIsDir),
                  PathVariable('BOOST_LIB', 'Boost library directory', BOOST_LIB, PathVariable.PathIsDir),
                 )

# Build environment setup
env = Environment(variables=vars)
env.Help("Spherical data partitioning and duplication utilities.")

if env.get('BOOST_INC'):
    BOOST_INC = env['BOOST_INC']
    if not os.path.isabs(BOOST_INC):
        BOOST_INC = os.path.normpath(os.path.join(GetLaunchDir(), BOOST_INC))
    env.Append(CPPPATH=[BOOST_INC])
if env.get('BOOST_LIB'):
    BOOST_LIB = env['BOOST_LIB']
    if not os.path.isabs(BOOST_LIB):
        BOOST_LIB = os.path.normpath(os.path.join(GetLaunchDir(), BOOST_LIB))
    env.Append(LIBPATH=[BOOST_LIB])
    env.Append(RPATH=[BOOST_LIB])
env.Append(CCFLAGS=['-g', '-Wall', '-Wextra'])
if env['debug']:
    variant_dir = 'build/debug'
else:
    variant_dir = 'build/release'
    env.Append(CCFLAGS=['-O2'])
if not os.path.isabs(env['prefix']):
    # `scons prefix=<dir> install` seems buggy when <dir> is relative
    env['prefix'] = os.path.normpath(os.path.join(GetLaunchDir(), env['prefix']))

# Check for required boost libraries
boostLibs = ['system',
             'thread',
             'filesystem',
             'program_options',
            ]

def checkBoost(lib, conf):
    return (conf.CheckLib('boost_' + lib + '-mt', language='C++') or
            conf.CheckLib('boost_' + lib, language='C++'))

if not GetOption('help'):
    if not GetOption('clean'):
        conf = Configure(env)
        if filter(lambda l: not checkBoost(l, conf), boostLibs):
            print >> sys.stderr, 'Missing required boost library!'
            Exit(1)
        env=conf.Finish()

env.Alias('install', ['$prefix/bin', '$prefix/ups'])

# Recurse
if not GetOption('help'):
    SConscript(dirs='.', variant_dir=variant_dir, exports='env')
