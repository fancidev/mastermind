#!/usr/bin/perl
#
# Run tests for the mmstrat utility.
#

use strict;
use warnings;

# Path to executable.
my $exec = 'mmstrat';
if ($#ARGV >= 0)
{
	$exec = $ARGV[0];
}

# Run the -v command to check that the executable exists.
print "$exec -v\n";
system($exec, "-v") == 0 or exit $?;

# Counter for test number.
my $number = 0;
my $failed = 0;

# Test strategies. Each test is specified by its command line arguments
# and the expected output. The "-S -q" switches are automatically appended
# to each entry. Individual options for specific strategies are also tested.
my @test_cases = (

	# Test heuristic strategies for standard Mastermind rules.
	"-r mm -s simple",          "6962:8:4",
	"-r mm -s minmax",          "5778:5:663",
	"-r mm -s minmax -no",      "5778:5:663",
	"-r mm -s minmax -nc",      "5780:5:663",
	"-r mm -s minmax -po",      "5797:6:44",
	"-r mm -s minavg",          "5696:6:3",
	"-r mm -s minavg -no",      "5696:6:3",
	"-r mm -s minavg -nc",      "5702:6:3",
	"-r mm -s minavg -po",      "5722:6:43",
	"-r mm -s entropy",         "5719:6:18",
	"-r mm -s entropy -no",     "5719:6:18",
	"-r mm -s entropy -nc",     "5726:6:12",
	"-r mm -s entropy -po",     "5787:6:58",
	"-r mm -s parts",           "5668:6:7",
	"-r mm -s parts -no",       "5668:6:7",
	"-r mm -s parts -nc",       "5684:6:7",
	"-r mm -s parts -po",       "5714:7:2",

	# Test optimal strategies.
	"-r mm -s optimal",         "5625:6:7",
	"-r mm -s optimal -O 1",    "5625:6:7",
	"-r mm -s optimal -po",     "5629:6:7",
	"-r bc -s optimal -po",     "26374:7:126",

	# Test -md switch for optimal strategies.
	#"-r mm -s optimal -md 10",  "5625:6:7",
	#"-r mm -s optimal -md 6",   "5625:6:7",
	#"-r mm -s optimal -md 5",   "5626:5:x",
	#"-r mm -s optimal -md 4",   "0:0:0",

	# Test different equivalence filters.
	"-r mm -s minavg -e default",    "5696:6:3",
	"-r mm -s minavg -e constraint", "5696:6:3",
	"-r mm -s minavg -e color",      "5696:6:3",
	"-r mm -s minavg -e none",       "5696:6:3",

	# Build strategy using 2 threads.
	"-r mm -mt 2 -s minmax",    "5778:5:663",
	"-r mm -mt 2 -s minavg",    "5696:6:3",
	"-r mm -mt 2 -s entropy",   "5719:6:18",
	"-r mm -mt 2 -s parts",     "5668:6:7",
	"-r mm -mt 2 -s optimal",   "5625:6:7",

	# Test Bulls and Cows rule for selected strategies.
	"-r bc -s simple",          "27511:8:41",
	"-r bc -s minmax",          "27030:7:181",
	"-r bc -s minavg",          "26551:7:87",
	"-r bc -s entropy",         "26409:8:1",
	"-r bc -s parts",           "26751:8:3",
	
	# Test other rules.
	"-r lg -mt 2 -s minavg",    "180287:7:789",
	"-r p1c2r -s optimal",      "3:2:1",
	"-r p3c9r -s optimal",      "3596:7:3",
	"-r p4c8n -s entropy",      "7880:7:2",
);

my $last_is_ok = 1;
for (my $i = 0; $i < $#test_cases; $i += 2)
{
	my $opt = $test_cases[$i];
	my $expect = $test_cases[$i+1];

	#print "Running test ", sprintf("%2d", ++$number), " ... ";
	print sprintf("Running test [ %3d / %3d ] ... ", ++$number, ($#test_cases+1)/2);

	my $args = "-S -q $opt";
	my $cmd = "$exec $args";
	my $actual = `$cmd`;
	chomp($actual);

	if ($actual eq $expect)
	{
		#print "OK\n";
		print "\r";
		$last_is_ok = 1;
	}
	else
	{
		print "FAILED\n";
		print "    Args:   $args\n";
		print "    Expect: $expect\n";
		print "    Actual: $actual\n";
		++$failed;
		$last_is_ok = 0;
	}
}

# Display summary.
print "\n" if $last_is_ok;
if ($failed == 0)
{
	print "All tests passed.\n";
	exit 0;
}
else
{
	print "$failed tests failed out of $number.\n";
	exit 1;
}
