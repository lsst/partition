# -*- python -*-
import os, sys
import SCons.Script

env = Environment()
opts = SCons.Script.Variables("custom.py")
opts.AddVariables(
          (PathVariable('prefix', 'Installation prefix', "build", PathVariable.PathIsDirCreate))
)

opts.Update(env)

# TODO: use Variables instead?
AddOption('--debug-build', dest='debug_build', action='store_true',
          default=False, help='Debug build')
AddOption('--boost-includes', dest='boost_includes', metavar='DIR',
          help='Boost include directory')
AddOption('--boost-libs', dest='boost_libs', metavar='DIR',
          help='Boost library directory')

env['ENV']['PATH'] = os.environ['PATH']

env.Append(CPPPATH=[GetOption('boost_includes')])
env.Append(LIBPATH=[GetOption('boost_libs')])
env.Append(RPATH=[GetOption('boost_libs')])
env.Append(CCFLAGS=['-g', '-Wall', '-Wextra'])
if GetOption('debug_build'):
    variant_dir = 'build/debug'
else:
    variant_dir = 'build/release'
    env.Append(CCFLAGS=['-O2'])
env.Alias('install', ['$prefix/bin', '$prefix/ups'])
env.Install('$prefix/ups', 'ups/partition.table')

env.Help("Spherical data partitioning and duplication utilities.")

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
    SConscript(dirs='.', variant_dir=variant_dir, exports='env')