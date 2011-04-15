use 5.008008;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME              => 'Math::Assistant',
    VERSION_FROM      => 'lib/Math/Assistant.pm',
    PREREQ_PM         => {
            'Algorithm::Combinatorics'   => '0.26',
            'Math::BigInt'	=> '1.993',
            'Math::BigFloat'	=> '1.993',
            'Carp'		=> '0',
            'Test::More'	=> '0',
    },
    dist => {
        COMPRESS => 'gzip -9f',
        SUFFIX => 'gz',
    },
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Math/Assistant.pm',
       AUTHOR         => 'Alessandro Gorohovski <angel@feht.dgtu.donetsk.ua>') : ()),
);
