#!/usr/bin/perl -w
# Computing the Number of Paths in a Grid Graph
# Hiroaki Iwashita <iwashita@erato.ist.hokudai.ac.jp>
# Copyright (c) 2013 ERATO MINATO Project
# $Id$

use strict;
use bigint;

my @cmd = qw(ggcount_multi ggcount_single);

my $dir = $0;
$dir =~ s/[^\/]*$//;
my $GGCOUNT;

for my $cmd (@cmd) {
    $GGCOUNT = $dir . $cmd;
    last if -X $GGCOUNT;
}

-X $GGCOUNT or die "Can't find " . join(" or ", @cmd) . "!\n";

my $USAGE = <<";;;";
Usage: $0 [<option>...] <size>
Options
  -64 : Use 64bit integer (default)
  -32 : Use 32bit integer
  -16 : Use 16bit integer
  -8  : Use 8bit integer
  -c  : Count cycles instead of paths
  -h  : Count Hamiltonian paths/cycles
  -v  : Print verbose messages
  -q  : Work quietly
;;;

$| = 1;

my $size       = -1;
my $maxModulus = 2**64;
my @opt;
my @arg;

for (@ARGV) {
    if (/^-[chvq]$/) {
        push(@opt, $_);
    }
    elsif (/^-(64|32|16|8)$/) {
        $maxModulus = 2**$1;
    }
    elsif (/^\d+$/ and $size < 0) {
        $size = $_ + 1;
    }
    else {
        push(@arg, $_);
    }
}

$size >= 2 and @arg == 0 or die $USAGE;

my @M;    # array of coprimes
my @A;    # array of [result, modulus]

sub gcd {
    my ($m, $n) = @_;
    return $n if $n == $m;
    ($n, $m) = ($m, $n) if $n < $m;

    while () {
        $n %= $m;
        return $m if $n == 0;
        return 1  if $n == 1;
        ($n, $m) = ($m, $n);
    }
}

sub exgcd {    # am + bn = gcd(m, n)
    my ($r0, $r1) = @_;
    my ($a0, $a1, $b0, $b1) = (1, 0, 0, 1);

    while ($r1 > 0) {
        my $q = int($r0 / $r1);
        my $r = $r0 % int($r1);
        my $a = $a0 - $q * $a1;
        my $b = $b0 - $q * $b1;
        $r0 = $r1, $r1 = $r;
        $a0 = $a1, $a1 = $a;
        $b0 = $b1, $b1 = $b;
    }

    ($a0, $b0, $r0);
}

sub isCoprime {
    my ($k) = @_;
    for my $m (@M) {
        return 0 if gcd($k, $m) != 1;
    }
    1;
}

sub getModulus {
    my ($i) = @_;
    my $k0 = @M ? $M[-1] - 1 : $maxModulus;

    for (my $k = $k0; $#M < $i; --$k) {
        $k > 1 or die "@M ... No more coprime number!\n";
        push(@M, $k) if isCoprime($k);
    }

    return $M[$i];
}

sub getAnswer {
    my $ans = 0;
    my $mmm = 1;

    for my $A (@A) {
        my ($n, $m) = @$A;
        $mmm *= $m;
    }

    for my $A (@A) {
        my ($n, $m) = @$A;
        my $mm = $mmm / $m;

        my ($a, $b, $c) = exgcd($mm, $m);
        $c == 1 or die "??? exgcd($mm, $m) = ($a, $b, $c)";

        $ans += $a * $mm * $n;
    }

    $ans %= $mmm;

    for my $A (@A) {
        my ($n, $m) = @$A;
        $ans % $m == $n or die "???";
    }

    ($ans, $mmm);
}

sub readLog {
    my ($logfile) = @_;
    open(my $log, "<", $logfile) or return undef;
    local $/;
    local $_ = <$log>;
    /^(\d+) *\(mod ((0x[\da-fA-F]+)|\d+)\)$/m or return undef;
    ($1, $3 ? hex($3) : $2);
}

my $upperBound = int do {
    no bigint;
    1.7817**($size * $size);
};

my $numRun = 0;
{
    printf STDERR "***** G%dx%d *****\n", $size, $size;
    my $mmm = 1;
    while ($mmm <= $upperBound) {
        $mmm *= getModulus($numRun++);
    }
    printf STDERR "Num of runs: %d\n",    $numRun;
    printf STDERR "Upper bound: %.10g\n", $upperBound;
}

my $answer;

for (my $i = 0;; ++$i) {
    my $m = getModulus($i);
    my @cmd = ($GGCOUNT, @opt, $size - 1, "%" . $m);

    my $name = join("_", @cmd);
    $name =~ s/^.*\///;
    my $log     = "$name.log";
    my $logging = "$name.logging";

    my ($n, $mm) = readLog($log);

    if (defined $n) {
        printf STDERR "  Log #%d: %s\n", $i + 1, join(" ", @cmd);
    }
    elsif (-f $logging) {
        printf STDERR "  Run #%d: %s\n", $i + 1, join(" ", @cmd);
        printf STDERR "  -> logging\n";
        next;
    }
    else {
        last if $i > $numRun;

        printf STDERR "  Run #%d: %s\n", $i + 1, join(" ", @cmd);
        my $pid = fork;
        defined $pid or die "fork failed: $!";

        if ($pid == 0) {
            rename($logging, "$name.faillog");
            open(STDOUT, ">", $logging) or die "$logging: Can't open for write";
            open(STDERR, ">&STDOUT") or die "Can't dup STDOUT: $!";
            exec(@cmd);
            exit;    # NOTREACHED
        }

        wait == $pid or die "lost the child!";
        $? == 0 or do {
            rename($logging, "$name.faillog");
            die "Exited abnormally\n";
        };

        rename($logging, $log);
        ($n, $mm) = readLog($log);
    }

    defined $n or do {
        print STDERR "$log: Can't get the result\n";
        next;
    };

    $mm == $m or die "Unexpected modulus in the log";
    push(@A, [$n, $m]);

    my $prev_a = $answer;
    my $modulus;
    ($answer, $modulus) = getAnswer();
    printf STDERR "  -> %.10g (mod %.10g)\n", $answer, $modulus;
    print $answer, "\n" if $answer == $prev_a;
}
