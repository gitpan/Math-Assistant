#use Test::More tests => 5;
use Test::More 'no_plan';

BEGIN { use_ok('Math::Assistant') };

use Math::Assistant qw(:algebra);

my $M = [ [1,4,2,0,-3], [2,9,5,2,1], [1,3,1,-2,9], [3,12,6,0,-8], [2,10,6,4,7] ];
is(Rank($M),3,'Rank of matrix 5x5');

$M = [
	[1,1,1,1,1,1,1],
	[64,32,16,8,4,2,1],
	[729,243,81,27,9,3,1],
	[4096,1024,256,64,16,4,1],
	[5**6, 5**5, 5**4, 5**3, 5**2, 5, 1],
	[6**6, 6**5, 6**4, 6**3, 6**2, 6, 1],
	[7**6, 7**5, 7**4, 7**3, 7**2, 7, 1],
    ];

is(Rank($M),7,'Rank of matrix 7x7');

is(Det($M),-24883200,'Determinant of matrix 7x7');

my $B = [ 7,127,1093,5461,19531,55987,137257 ];

is( scalar(grep $_ eq '1', @{ Solve_Det($M, $B) } ),7,'Solve an equation system');

$M = [ [1.3, 2.5, 3.1, 4.2], [2.1, 4.7, 6.2, 8.7], [34, 8.2, 12, 16], [78,16,24,33] ];
$B = [ 1.1, 2.2, 3.3, 4.4 ];

is(Det($M),34.754,'Determinant of float matrix 4x4');

my $solve = [ '-20427/173770', '-47993/34754', '317493/86885', '-27401/17377' ];

my $i = 0;
for my $s ( @{ Solve_Det($M, $B) } ) {
    is($s, $solve->[$i++], 'X');
}

