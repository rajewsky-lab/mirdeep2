use strict;
use warnings;
use Module::Build;

use 5.008;

my $extdeps;
$extdeps = [
	{	"dep" => "bowtie",
		"test" => "command -v -- bowtie --help",
		"hint" => "run install.pl" 
	},{	"dep" => "ViennaRNA",
		"test" => "command -v -- RNAfold -h",
		"hint" => "run install.pl"
	},{	"dep" => "randfold",
		"test" => "command -v -- randfold", # the actual command exits with 1 with no working args
		"hint" => "run install.pl to install squid and randfold "
	}
];

for my $dep (@{$extdeps}){
	warn "testing for ".$dep -> {'dep'}.".\n";
	CmdRunner($dep -> {'test'});
}


sub CmdRunner {
        my $ret;
        my $cmd = join(" ",@_);

        warn localtime( time() ). " [INFO] system call:'". $cmd."'.\n";
	
	#adding a `set -e; set -o pipefail` although safer wil crash the travis build env. 
        @{$ret} = `( $cmd )2>&1`;
	my $exitcode = $?;
        if ($exitcode == -1) {
                die localtime( time() ). " [ERROR] failed to execute: $!\n";
        }
	elsif ($exitcode & 127) {
                die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
                 ($exitcode & 127),  ($exitcode & 128) ? 'with' : 'without';
        }
	elsif ($exitcode != 0) {
                die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
                 ($exitcode & 127),  ($exitcode & 128) ? 'with' : 'without';
        }else {
               	warn localtime( time() ). " [INFO] " . sprintf "child exited with value %d\n", $exitcode >> 8;
        }
	return @{$ret};
}


my @exe_files = glob("src/*.pl");
my $builder = Module::Build->new(
    module_name         => 'Bio::miRDeep2',
    license             => 'perl',
    dist_author         => 'Sebastian Mackowiak & Marc Friedländer',
    dist_abstract       => 'miRDeep2 discovers active known or novel miRNAs from deep sequencing data.',
    create_makefile_pl  => 0,     #'traditional',
    script_files        => [@exe_files],
    create_readme       => 1,
    requires            => {
        'PDF::API2'      => '>=2.028',
    },
    build_requires      => {
        'Test::More'          => '0.47',
    },
    meta_merge => {
        resources => {
            repository => 'https://github.com/rajewsky-lab/mirdeep2',
            bugtracker => 'https://github.com/rajewsky-lab/mirdeep2/issues'
        }
    },
);

$builder->create_build_script();
