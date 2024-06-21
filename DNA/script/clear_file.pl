#!/usr/bin/perl -w
use strict;
die "$0 files list" if @ARGV < 1 ;
unlink @ARGV;
