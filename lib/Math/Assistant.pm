package Math::Assistant;

use 5.008008;
use strict;
use warnings;

use Carp;
use Algorithm::Combinatorics qw( permutations );
use Math::BigInt qw( bgcd );

require Math::BigFloat;
require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [ qw(
			Rank Det Solve_Det Gaussian_elimination
			) ],
		    'algebra' => [ qw( Rank Det Solve_Det ) ],
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw( );

our $VERSION = '0.0417';

use base qw(Exporter);


# Rank of any integer matrix
# input:
#	$Dimention
# return:
#	$rank (order)
sub Rank{
    my $M = &Gaussian_elimination; # triangular matrix
    my $rank = 0;

    for( @$M ){ # rows
	for( @{$_} ){ # element in row
	    $rank++, last if $_
	}
    }
    $rank;
}


# Method of Gaussian elimination
# input:
#	$Dimention
# return:
#	$elimination_Dimention
sub Gaussian_elimination{
    my $M_ini = shift;

    my $t = &_test_matrix( $M_ini );
    if( $t > 3 ){
	croak("Use of uninitialized value in matrix");

    }elsif( $t > 1 ){
	croak("Bad matrix");
    }

    my $rows = $#{$M_ini}; # number of rows
    my $cols = $#{$M_ini->[0]}; # number of columns

    my %zr; # for rows with 0
    my %zc; # for columns with 0

    my $M; # copy
    # search of rows with 0
    for( my $i = 0; $i <= $rows; $i++ ){
	for( my $j = 0; $j <= $cols; $j++ ){
	    $M->[$i][$j] = $M_ini->[$i][$j];
	    unless($M->[$i][$j]){ # zero element
		$zr{$i}++;
		$zc{$j}++;
	    }
	}
    }

    # Check Float (Fraction)
    for my $i ( @$M ){
	my $max_frac = &_max_len_frac($i, 0);

	if($max_frac){
	    $_ *= 10**$max_frac for @$i;
	}
    }

    if(keys %zr){ # yes, zero rows (columns)

	# I supplement indexes of nonzero rows
	for( my $i = 0; $i <= $rows; $i++ ){
	    $zr{$i} += 0;
	}

	my $R; # temp copy of matrix M
	my $v = 0;
	# replacement of rows
	for my $i ( sort { $zr{$a} <=> $zr{$b} } keys %zr ){
	    for( my $j = 0; $j <= $cols; $j++ ){
		$R->[$v][$j] = $M->[$i][$j];
	    }
	    $v++;
	}

	# I supplement indexes of non zero columns
	for( my $j = 0; $j <= $cols; $j++ ){
	    $zc{$j} += 0;
	}

	$v = 0;
	# replacement of columns
	for my $j ( sort { $zc{$b} <=> $zc{$a} } keys %zc ){
	    for( my $i = 0; $i <= $rows; $i++ ){
		$M->[$i][$v] = $R->[$i][$j];
	    }
	    $v++;
	}
    }

    undef %zr;
    undef %zc;

    for(my $n = 0; $n < $rows; $n++){
	last if $n > $cols; # other zero rows

	# Replacement of zero diagonal element
        my $s = 0;
M_Gauss1:
	while( $M->[$n][$n] == 0 && $s < $rows - $n ){
	    $s++;
	    for( my $j = $n + 1; $j <= $cols; $j++ ){

		if( $M->[$n][$j] ){ # no zero element
		    # shift of columns $n <-> $j
		    for( my $i = 0; $i <= $rows; $i++ ){
			($M->[$i][$j], $M->[$i][$n]) = ($M->[$i][$n], $M->[$i][$j]);
		    }
		    last M_Gauss1;
		}
	    }

	    # all zero elements in row
	    for( my $j = $n; $j <= $cols; $j++ ){
		# shift of rows $n <-> $n+1
		($M->[$n][$j], $M->[$n+$s][$j]) = ($M->[$n+$s][$j], $M->[$n][$j]);
	    }
	}

	last unless $M->[$n][$n]; # zero elements of rows

	for( my $i = $n+1; $i <= $rows; $i++ ){

	    # Divisibility check totally
	    my($k,$d);
	    my $b = $M->[$i][$n] / $M->[$n][$n];
	    if($b == int($b)){
		$k = -$b;
		$d = 1;
	    }else{
		$k = -$M->[$i][$n];
		$d = $M->[$n][$n];
	    }

	    # column (element in row)
	    for( my $j = $n; $j <= $cols; $j++ ){
		$M->[$i][$j] = $M->[$i][$j]*$d + $k*$M->[$n][$j];
	    }
	}

	my $gcd = bgcd( @{$M->[$n]}[$n..$cols] );
	if($gcd > 1){
	    $_ = $_ / $gcd for @{$M->[$n]}[$n..$cols];
	}
    }
    $M;
}


# Interger determinant for quadratic matrix
# input:
#	$Dimention
# return:
#	determinant
#	undef
sub Det{
    my $M_ini = shift;

    my $t = &_test_matrix($M_ini);
    if( $t > 3 ){
	croak("Use of uninitialized value in matrix");

    }elsif( $t ){
	croak("Matrix is not quadratic");
    }

    my $det = 0;

    my $rows = $#{$M_ini}; # number of rows
    my $cols = $#{$M_ini->[0]}; # number of columns

    return $M_ini->[0][0] if $#{$M_ini} < 1; # dim matrix = 1

    # Check Float (Fraction)
    my $max_frac = 0;
    for my $i ( @$M_ini ){
	$max_frac = &_max_len_frac($i, $max_frac);
    }

    my $M = $max_frac ? [ map{ [ map{ $_ * 10**$max_frac } @$_ ] } @$M_ini ] : $M_ini;

    my $iter = permutations([(0..$#{$M})]);
    while(my $prm = $iter->next){
	my $sign = 1; # +/-
	my $p = 1;
	my $i = 0;
	my $j = 0;
	for my $x ( @$prm ){
	    $p *= $M->[$i++][$x];

	    # Inversions in permutations
	    for( my $k = ++$j; $k <= $#{$prm}; $k++ ){
		$sign *= -1 if $x > $prm->[$k];
	    }
	}
	$det += $sign * $p;
    }
    return $max_frac ? $det/10**($max_frac * ($rows + 1)) : $det;
}


sub Solve_Det{
    my $M = shift || croak("Missing matrix");
    my $B = shift || croak("Missing vector");
    my $opts = shift;

    my $rows = $#{$M}; # number of rows
    my $cols = $#{$M->[0]}; # number of columns

    my $t = &_test_matrix($M);
    if( $t > 3 ){
	croak("Use of uninitialized value in matrix");

    }elsif( $t ){
	croak("Matrix is not quadratic");
    }

    croak("Vector doesn't correspond to a matrix") if $rows != $#{$B};
    croak("Use of uninitialized value in vector") if scalar( grep ! defined $_, @$B );

    if(defined $opts){
	if(exists $opts->{'eqs'}){
	    die "Unknown parameter \'$opts->{'eqs'}\'!\n"
		unless $opts->{'eqs'}=~/^(?:row|column)/i;

	}else{
	    $opts->{'eqs'} = 'row';
	}
    }else{
	$opts->{'eqs'} = 'row';
    }

    my $solve;

    # main determinant
    my $det_main = &Det( $M ) || return undef; # no one solution

    for( my $v = 0; $v <= $cols; $v++ ){

	my $R; # copy of matrix M
	for( my $i = 0; $i <= $rows; $i++ ){

	    if($opts->{'eqs'}=~/^col/i){
		$R->[$i] = $v == $i ? $B : $M->[$i];

	    }else{
		for( my $j = 0; $j <= $cols; $j++ ){
		    $R->[$i][$j] = $v == $j ? $B->[$i] : $M->[$i][$j];
		}
	    }

	}

	my $det = &Det( $R );
	my $dm = $det_main;

	if( $det ){
	    if( $dm < 0 ){
		$dm *= -1;
		$det *= -1;
	    }

	    # Check Float (Fraction)
	    my $max_frac = &_max_len_frac([$dm, $det], 0);

	    if($max_frac){
		$_ *= 10**$max_frac for($dm, $det);
	    }

	    my $gcd = bgcd(abs($det), abs($dm));
	    if($gcd > 1){
		$det = $det / $gcd;
		$dm = $dm / $gcd;
	    }

	    $solve->[$v] = $dm == 1 ? $det : "$det/$dm";

	}else{
	    $solve->[$v] = 0;
	}
    }
    $solve;
}


sub _max_len_frac{
    my($M, $max_frac) = @_;

    for( @$M ){
	next if Math::BigInt::is_int($_);
	my(undef, $frac) = Math::BigFloat->new($_)->length();
	$max_frac = $frac if $frac > $max_frac;
    }
    $max_frac;
}


sub _test_matrix{
    my $M = shift;

    return 4 if scalar( grep{ grep ! defined $_, @{$_} } @$M );

    my $rows = $#{$M}; # number of rows
    my $cols = $#{$M->[0]}; # number of columns

    my $quadra = 0;
    $quadra = scalar( grep $#{$_} != $rows, @$M );
    my $reqtan = 0;
    $reqtan = scalar( grep $#{$_} != $cols, @$M );

    my $res = 0;
    $res++ if $quadra;
    $res += 2 if $reqtan;

    $res;
}

1;
__END__

=head1 NAME

Math::Assistant - functions for various exact algebraic calculations

=head1 SYNOPSIS

  use Math::Assistant qw(:algebra);

  my $M = [ [4,1,4,3], [3,-4,7,5], [4,-9,8,5], [-3,2,-5,3], [2,2,-1,0] ];

  # Rank of rectangular matrix
  my $rank = Rank($M);
  print "Rank = $rank\n";

  shift @{$M}; # Now a quadratic matrix

  # Determinant of quadratic matrix
  my $determinant = Det($M);
  print "Determinant = $determinant\n";

  # Solve an equation system
  my $B = [ 1,2,3,4 ];
  my $solve = Solve_Det($M, $B); # eqs => 'row' (default)
  print "Equations is rows of matrix:\n";
  print "$_\n" for @{$solve};

  use Math::BigRat;
  print(Math::BigRat->new("$_")->numify(),"\n") for @{$solve};

  print "Equations is columns of matrix:\n";
  print "$_\n" for @{ Solve_Det( $M, $B, {'eqs' => 'column' } ) ;


will print

    Rank = 4
    Determinant = -558
    Equations is rows of matrix:
	433/279
	-32/279
	-314/279
	70/93

    1.55197132616487
    -0.114695340501792
    -1.12544802867384
    0.752688172043011

    Equations is columns of matrix:
	283/186
	-77/93
	11/62
	13/93


=head1 DESCRIPTION

The module contains important algebraic operations: matrix rank, determinant and
solve an equation system. The integer values are accepted.
Calculations with the raised accuracy.


=head1 SUBROUTINES

Math::Assistant provides these subroutines:

    Rank(\@matrix)
    Det(\@matrix)
    Solve_Det(\@A_matrix, \@b_vector, {eqs => 'row|column' } )
    Gaussian_elimination(\@matrix)

All of them are context-sensitive.


=head2 Rank(\@matrix)

Calculates rank of rectangular (quadratic or non-quadratic) C<@matrix>.
Rank is a measure of the number of linear independent row and column
(or number of linear independent equations in the case of a matrix representing
an equation system).


=head2 Det(\@matrix)

This subroutine returns the determinant of the C<@matrix>.
Only quadratic matrices have determinant. For non-quadratic matrix returns C<undef>.


=head2 Solve_Det(\@A_matrix, \@b_vector, {eqs => 'row|column' } )

Use this subroutine to actually solve an equation system.

Matrix "C<@A_matrix>" must be quadratic matrix of your equation system
C<A * x = b>. By default the equations are in rows.

The input vector "C<@b_vector>" is the vector "b" in your equation system
C<A * x = b>, which must be a row vector and have the same number of
elements as the input matrix "C<@A_matrix>" have rows (columns).

The subroutine returns the solution vector "C<$x>"
(which is the vector "x" of your equation system C<A * x = b>) or C<undef>
is no solution.

    # Equation system:
    # x1 + x2 + x3 + x4 + x5 + x6 + x7 = 4
    # 64*x1 + 32*x2 + 16*x3 + 8*x4 + 4*x5 + 2*x6 + x7 = 85
    # ...................................
    # 7**6*x1 + 7**5*x2 + 7**4*x3 + 7**3*x4 + 7**2*x5 + 7*x6 + x7 = 120100

    my $M = [
	[1,1,1,1,1,1,1],
	[64,32,16,8,4,2,1],
	[729,243,81,27,9,3,1],
	[4**6, 4**5, 256,64,16,4,1],
	[5**6, 5**5, 5**4, 5**3, 5**2, 5, 1],
	[6**6, 6**5, 6**4, 6**3, 6**2, 6, 1],
	[7**6, 7**5, 7**4, 7**3, 7**2, 7, 1],
	 ];

    my $B = [ 4, 85, 820, 4369, 16276, 47989, 120100 ];

    print "$_\t" for @{ Solve_Det($M, $B) };

will print:

    1 0 1 0 1 0 1

Other example:

    # Equation system:
    # 1.3*x1 + 2.1*x2 + 34*x3 + 78*x4 = 1.1
    # 2.5*x1 + 4.7*x2 + 8.2*x3 + 16*x4 = 2.2
    # 3.1*x1 + 6.2*x2 + 12*x3 + 24*x4 = 3.3
    # 4.2*x1 + 8.7*x2 + 16*x3 + 33*x4 = 4.4

    $M = [  [1.3, 2.5, 3.1, 4.2],
	    [2.1, 4.7, 6.2, 8.7],
	    [34,  8.2, 12,  16],
	    [78,  16,  24,  33] ];
    print "$_\t" for @{ Solve_Det($M, [ 1.1, 2.2, 3.3, 4.4 ], {eqs => 'column'} ) };

will print:

    -38049/17377  22902/17377  35101/34754  -36938/86885


=head2 Gaussian_elimination(\@matrix)

This subroutine returns matrix Gaussian elimination of the C<@matrix>.
The initial C<@matrix> does not vary.


=head1 EXPORT

Math::Assistant exports nothing by default.
Each of the subroutines can be exported on demand, as in

  use Math::Assistant qw( Rank );

the tag C<algebra> exports the subroutines C<Rank>, C<Det>, C<Solve_Det>:

  use Math::Assistant qw(:algebra);

and the tag C<all> exports them all:

  use Math::Assistant qw(:all);


=head1 DEPENDENCIES

Math::Assistant is known to run under perl 5.8.8 on Linux.
The distribution uses L<Algorithm::Combinatorics> for the subroutine C<Det()>,
L<Math::BigInt>, L<Math::BigFloat> and L<Carp>.


=head1 SEE ALSO

L<Math::MatrixReal> is a Perl module that offers similar features.


=head1 AUTHOR

Alessandro Gorohovski, E<lt>angel@domashka.kiev.uaE<gt>, E<lt>angel@feht.dgtu.donetsk.uaE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2012 by A. N. Gorohovski

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.


=cut
