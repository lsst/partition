# -*- python -*-
Import('env')

SConscript(dirs=['src', 'tests'], exports='env')
if env.Dir('$prefix/ups').abspath != env.Dir('#ups').abspath:
    env.Install('$prefix/ups', '#ups/partition.table')
