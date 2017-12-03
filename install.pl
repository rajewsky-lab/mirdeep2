#!/usr/bin/env perl

use strict;
use warnings;

print STDERR "\n
###############################################################################################
# Last update: Dec 03, 2017
# This is the miRDeep2 installer
# It is tested under a bash and zsh shell
# It will try to download all necessary third-party tools and install them. 
# To run this installer you need to have mirdeep2.0.0.8 or later
#
###############################################################################################

";
my $dir=`pwd 2>&1`;
chomp $dir;

if(not $ENV{PERL5LIB}){$ENV{PERL5LIB}="$dir/lib/perl5";}
if($ENV{PERL5LIB} !~ /$dir\/lib\/perl5/){
    $ENV{PERL5LIB}.=":$dir/lib/perl5";
}

my $time =time;
my $new='no';

$new = $ARGV[0] if($ARGV[0]);
my $shell=$ENV{'SHELL'};
my $shellconf='.bash_profile';
if($shell =~ /zsh/){
	if(-f "$ENV{'HOME'}/.zshenv"){
		$shellconf='.zshenv';
	}elsif(-f "$ENV{'HOME'}/.zshrc"){
		$shellconf='.zshrc';
	}else{
		die "Could not determine file for setting environment variables\n";
	}
}



if($new !~ /no/){
	print STDERR "making backup of .bashrc,$shellconf and .cshrc and removing all entries of mirdeep in those files\n";
	rem_mirdeep(".bashrc");
	rem_mirdeep("$shellconf");
	rem_mirdeep(".cshrc");
}

my $grep;
$grep=`which grep`;
if(not $grep){
	$grep=`which ggrep`;
}
if(not $grep){
	die "No grep found on system\n";
}
chomp $grep;


my $gcc=`gcc --version 2>&1`;
if($gcc !~ /(GCC)/i and $gcc !~ /clang/i){
	print STDERR "\nError:\n\tno gcc compiler installed. Please install a gcc compiler\n";
	my $r=`uname -s`;
	chomp $r;
	if($r !~ /Linux/i){
		print STDERR "==>     If you are using MacOS then you probably need to install Xcode with commandline-tools from the appstore!\n\n";
	}
	exit;
}else{
	if($gcc =~ /^gcc\s*\S*\s*(\d+\S+)\s*/){
		print STDERR "gcc version: $1                                      already installed, nothing to do ...\n";
	}
	if($gcc =~ /clang/){
		print STDERR "clang installed                                      already installed, nothing to do ...\n";
	}
} 

my %progs;
$progs{bowtie}=0;
$progs{RNAfold}=0;
$progs{randfold}=0;
$progs{zlib}=0;
$progs{pdf}=0;
$progs{ttf}=0;

my $wget=`wget 2>&1`;
my $curl=`which curl`;

my $dtool='';
my $dopt='';

if($wget =~ /URL/i){
	$dtool ="wget";

}elsif($curl){
	$dtool ="curl -L"; ## forces curl to follow redirections
	$dopt=" -O";
}else{
	die "No commandline download tool found on your system. Please install wget or curl on your machine\n";
}

if(not -d 'bin'){
	#creating bin directory which will also contain other executables in the end\n";
	my $ret=system("cp -r src bin");
	if(not $ret){ 
		print STDERR "bin directory created successful\n";
	}else{
		die "Could not create binary directory\n";
	}
}


my $err;
my $dfile='';

##only attach to config file if not yet existing
my $in=`$grep "$dir/bin" ~/.bashrc`;



## set install dir
my $install_bin_dir="$dir/bin";
foreach my $e(@ARGV){
	if($e =~ /install-dir=(.+)/){
		$install_bin_dir=$1;
	}
}

if($install_bin_dir ne "$dir/bin"){
	## check if it is writable and existent
	if(not -d $install_bin_dir){
		print STDERR "The given installation directory by argument install-dir is not existent\nexecutable files will be put into $dir/bin instead\n";
		$install_bin_dir="$dir/bin";
	}else{	
		chdir $install_bin_dir;
		my $ret=system("touch mirdeep_test_file");
		if(not $ret){
			system("rm mirdeep_test_file");
		}else{
			print STDERR "The given installation directory by argument install-dir is either not existent or not writeable\nexectable files will be put into $dir/bin instead\n";
			$install_bin_dir="$dir/bin";
		}
		chdir "$dir";
	}
}

## check if we have the install path in our files 
$in=`$grep "$install_bin_dir" ~/$shellconf`;
if(not $in){
	my $ret=`$grep $install_bin_dir ~/$shellconf |$grep PATH`;
	if(not $ret){
		`echo 'export PATH=\$PATH:$install_bin_dir' >> ~/$shellconf`;
	}
}

## add this temporarily to make perl installation possible on some systems
print STDERR "Checking environment variables ...\n";
my $g=`$grep PERL_MB_OPT ~/$shellconf \|$grep install_base `;
if($g){
}else{
	print STDERR "adding variables PERL_MB_OPT,PERL_MM_OPT,PERL5LIB to $shellconf\n";
	`echo >> ~/$shellconf`;
	`echo 'PERL_MB_OPT=\"--install_base $ENV{'HOME'}/perl5\";export PERL_MB_OPT' >> ~/$shellconf`;
	`echo 'PERL_MM_OPT=\"INSTALL_BASE=$ENV{'HOME'}/perl5\";export PERL_MM_OPT' >> ~/$shellconf`;
	$g=`grep $dir/lib/perl5 ~/$shellconf`;
	if(not $g){
		$g=`grep PERL5LIB ~/$shellconf`;
		if(not $g){
			`echo 'export PERL5LIB=$dir/lib/perl5' >> ~/$shellconf`;
		}else{
			`echo 'export PERL5LIB=\$PERL5LIB:$dir/lib/perl5' >> ~/$shellconf`;
		}	
	}

	`echo >> ~/$shellconf`;
	print STDERR "please run the install.pl script again in a new terminal window or just type

	source ~/$shellconf
	perl install.pl

	so that the new environment variables are visible to the install.pl script\n";

	exit;
}

$g=`$grep $dir/lib/perl5 ~/$shellconf`;
if(not $g){
	$g=`grep PERL5LIB ~/$shellconf`;
	if(not $g){
		`echo 'export PERL5LIB=$dir/lib/perl5' >> ~/$shellconf`;
	}else{
		`echo 'export PERL5LIB=\$PERL5LIB:$dir/lib/perl5' >> ~/$shellconf`;
	}	

	`echo >> ~/$shellconf`;
	print STDERR "please run the install.pl script again in a new terminal window or just type

	source ~/$shellconf
	perl install.pl

	so that the new environment variables are visible to the install.pl script\n";
	exit;
}

my $in2;
if(-f "~/.cshrc"){
	$in2=`$grep "$install_bin_dir" ~/.cshrc`;
	if(not $in2){
		`echo 'setenv PATH \$PATH:$install_bin_dir' >> ~/.cshrc`;
	}
}

my $binnew=1;

if(not -d "essentials"){
	`mkdir essentials`;
}else{
	$binnew=0;
}

chdir("essentials");

my $a=`uname -a`;

my $bowtie;

my $bowtie_version="1.1.1";
if($dtool =~ /curl/){
	`$dtool https://sourceforge.net/projects/bowtie-bio/files/bowtie/ > to_del`;
}else{
	`$dtool https://sourceforge.net/projects/bowtie-bio/files/bowtie/ -O to_del`;
}
open IN,"to_del" or die "No bowtie_file_info file found\n";
while(<IN>){
	if(/projects\/bowtie-bio\/files\/bowtie\/(\d\.\d+\.*\d*)\//){
		$bowtie_version=$1;
		last;
	}
}
close IN;

$bowtie_version="1.1.1"; #we force this version, cause 1.2 needs the TBB lib installed 

my $bv=$bowtie_version;
if($bv =~ /(\d\.\d)\.0/){
	$bv=$1;
}


my $ret=checkBIN("bowtie","Usage");
if($ret == 0){
	print STDERR "bowtie                                           already installed, nothing to do ...\n";
	$progs{bowtie} = 1;
}else{
	if(not -d "bowtie-$bowtie_version"){
		## this needed to be added cause the authors removed the 0 in the version number for the filename


		print STDERR "Downloading bowtie $bowtie_version binaries\n\n";
		if($a =~ /Darwin/i){ ## download mac version
			$bowtie = "bowtie-$bv-macos-x86_64.zip";
		}elsif($a =~ /x86_64/i){
			$bowtie = "bowtie-$bv-linux-x86_64.zip";
		}else{
			$bowtie = "bowtie-$bv-src.zip";
		}

		if(not -f $bowtie){
			if(check("http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/$bowtie_version/$bowtie")){
				$err=system("$dtool http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/$bowtie_version/$bowtie $dopt");

				if($err){
					die "\nError:\n\t$bowtie could not be downloaded\n\n\n";
				}
			}elsif(check("http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/old/$bowtie_version/$bowtie")){
				$err=system("$dtool http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/old/$bowtie_version/$bowtie $dopt");
				if($err){
					die "\nError:\n\t$bowtie could not be downloaded\n\n\n";
				}
			}elsif(check("https://sourceforge.net/projects/bowtie-bio/files/bowtie/$bowtie_version/$bowtie",$bowtie)){
				if(not -f $bowtie){
					$err=system("$dtool https://sourceforge.net/projects/bowtie-bio/files/bowtie/$bowtie_version/$bowtie $dopt");
				}else{
					$err=0;
				}
				if($err){
					die "\nError:\n\t$bowtie could not be downloaded\n\n\n";
				}

			}else{
				die "\nError:\n\t$bowtie not found on server http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/ \n\n\n";
			}
		}

		if(not -f "$bowtie"){
			die "$bowtie download failed \nPlease try to download bowtie manually from here http://bowtie-bio.sourceforge.net/index.shtml";
		}

		print STDERR "Installing bowtie binaries\n\n";
		$err=system("unzip -u $bowtie 1>> install.log 2>>install_error.log");

		if($err){
			die "unzip $bowtie was not successful\n";
		}
	}

	chdir "$install_bin_dir";

	buildgood("$dir/essentials/bowtie-$bv/bowtie","bowtie");

	if(not -f "bowtie"){
		system("ln -s $dir/essentials/bowtie-$bv/bowtie* .");
	}
	chdir "$dir/essentials/";

}

$ret = checkBIN("RNAfold -h","usage");

if($ret == 0){
	print STDERR "RNAfold                                          already installed, nothing to do ...\n";
	$progs{RNAfold}=1;
}else{
	if(not -d "ViennaRNA-1.8.4"){
		$dfile="ViennaRNA-1.8.4.tar.gz";
		if(not -f $dfile){
			print STDERR "Downloading Vienna package now\n\n";
			my $ptc="https://www.tbi.univie.ac.at/RNA/download/sourcecode/1_8_x/ViennaRNA-1.8.4.tar.gz";
			# if(check("http://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-1.8.4.tar.gz")){ # address changed 
			if(check($ptc,$dfile)){
				if(not -f $dfile){
					$err=system("$dtool $ptc $dopt");
				}else{
					$err=0;
				}
				if($err){
					die "Download of Vienna package not successful\n\n";
				}
			}else{
				die "Vienna package not found at  http://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-1.8.4.tar.gz
				Please try to download the Vienna package from here http://www.tbi.univie.ac.at/RNA/RNAfold.html 
				\n";
			}
		}

		if(not -f "ViennaRNA-1.8.4.tar.gz"){
			die "Vienna package download failed\n";
		}
	}

	if(not -f "ViennaRNA-1.8.4/Progs/RNAfold"){
		print STDERR "Installing Vienna package now \n\n";
		`tar xzf ViennaRNA-1.8.4.tar.gz`;
		chdir("ViennaRNA-1.8.4/");
		chdir("lib");

		open IN,"<fold.c" or die "File fold.c not found\n";
		open OUT,">fold.c.new" or die "Cannot generate file fold.c.new\n";
		while(<IN>){
			if(/inline\s+(int\s+LoopEnergy.+$)/i){
				print OUT "$1";
			}elsif(/^inline\s+(int\s+HairpinE.+$)/i){
				print OUT "$1";
			}else{
				print OUT;
			}
		}
		close OUT;

		`mv fold.c fold.c.orig`;
		`mv fold.c.new fold.c`;
		chdir("..");

		`./configure --prefix=$dir/essentials/ViennaRNA-1.8.4/install_dir`;
		`make 1>> ../install.log 2>> ../install_error.log`;
		`make install 1>> ../install.log 2>> ../install_error.log`;

		buildgood("$dir/essentials/ViennaRNA-1.8.4/install_dir/bin/RNAfold","RNAfold");

		chdir("..");
		chdir "$install_bin_dir";
		if(not -f "RNAfold"){
			system("ln -s $dir/essentials/ViennaRNA-1.8.4/install_dir/bin/RNAfold .");
		}
		chdir "$dir/essentials/"	
	}
}

#$in = `$grep "$dir/essentials/ViennaRNA-1.8.4/install_dir/bin:*" ~/.bashrc`;
#if(not $in){
#    print STDERR "Vienna package path has been added to \$PATH variable\n"; 
#    `echo 'export PATH=\$PATH:$dir/essentials/ViennaRNA-1.8.4/install_dir/bin' >> ~/.bashrc`;
#}

#$in = `$grep "$dir/essentials/ViennaRNA-1.8.4/install_dir/bin:*" ~/$shellconf`;
#if(not $in){
#    print STDERR "Vienna package path has been added to \$PATH variable\n"; 
#    `echo 'export PATH=\$PATH:$dir/essentials/ViennaRNA-1.8.4/install_dir/bin' >> ~/$shellconf`;
#}



#$in2 = `$grep "$dir/essentials/ViennaRNA-1.8.4/install_dir/bin:*" ~/.cshrc`;
#if(not $in2){
#`echo 'setenv PATH \$PATH:$dir/essentials/ViennaRNA-1.8.4/install_dir/bin' >> ~/.cshrc`;
#}




$ret = checkBIN("randfold","let7");

#my $randf = `randfold -h`;

#if($randf =~ /no\s*randfold/i){ ## this should work
if($ret == 0){
	print STDERR "randfold\t\t\t\t\t already installed, nothing to do ...\n";
	$progs{randfold}=1;
}else{

	$dfile="squid-1.9g.tar.gz";
	if(not -f $dfile){
		print STDERR "Downloading SQUID library now\n\n";
		`$dtool http://eddylab.org/software/squid/squid.tar.gz $dopt`;
		`mv squid.tar.gz $dfile`;
	}
	if(not -f $dfile){
		system("cp ../squid-1.9g.tar.gz squid-1.9g.tar.gz");
	}

	if(not -f "squid-1.9g.tar.gz"){
		die "squid could not be downloaded\n Please try to download the library from here http://selab.janelia.org/software.html";
	}

	if(not -d "squid-1.9g" and -f "squid-1.9g.tar.gz"){
		print STDERR "Extracting squid and configuring it now\n\n"; 
		`mkdir squid-1.9g`;
		`tar xzf squid-1.9g.tar.gz -C squid-1.9g`;
		my $a=`ls squid-1.9g`;
		chomp $a;
		`mv squid-1.9g/$a/* squid-1.9g/`;
		`rm -rf squid-1.9g/$a`;

		chdir("squid-1.9g");
		`./configure 1>>../install.log 2>>../install_error.log`;
		`make 1>>../install.log 2>>install_error.log`;

		buildgood("$dir/essentials/squid-1.9g/libsquid.a");

		chdir("..");
	}

	$dfile="randfold-2.0.tar.gz";
	if(not -f $dfile ){
		print STDERR "Downloading randfold now\n\n";
		`$dtool http://bioinformatics.psb.ugent.be/supplementary_data/erbon/nov2003/downloads/randfold-2.0.tar.gz $dopt`;
	}

	if(not -f "randfold-2.0.tar.gz"){
		die "randfold could not be downloaded\nPlease try to download randfold from here http://bioinformatics.psb.ugent.be/software/details/Randfold\n";
	}

	if(not -d "randfold-2.0" and -f "randfold-2.0.tar.gz"){
		print STDERR "Installing randfold now\n\n";
		`tar xzf randfold-2.0.tar.gz`;

		chdir("randfold-2.0");

		open IN,"<Makefile" or die "File Makefile not found\n\n";
		open OUT,">Makefile_new" or die "Makefile_new could not be created\n\n";

		while(<IN>){
			if(/INCLUDE=-I\.\s*/i){
				print OUT "INCLUDE=-I. -I$dir/essentials/squid-1.9g -L$dir/essentials/squid-1.9g/\n";
			}else{
				print OUT;
			}
		}
		close IN;
		close OUT;

		## added so we can make it run on MacOSX as well.

		open IN,"<fold.c" or die "File fold.c not found\n";
		open OUT,">fold.c.new" or die "Cannot generate file fold.c.new\n";
		while(<IN>){
			if(/^inline\s+(int\s+LoopEnergy.+$)/i){
				print OUT "$1";
			}elsif(/^inline\s+(int\s+HairpinE.+$)/i){
				print OUT "$1";
			}else{
				print OUT;
			}
		}
		close OUT;


		`mv fold.c fold.c.orig`;
		`mv fold.c.new fold.c`;

		`mv Makefile Makefile.orig`;
		`mv Makefile_new Makefile`;

		`make 1>>../install.log 2>>../install_error.log`;
		buildgood("$dir/essentials/randfold-2.0/randfold","randfold");
		chdir("..");
	}

#    $in = `$grep "$dir/essentials/randfold-2.0:*" ~/.bashrc`;
#    if(not $in){ 
#        print STDERR "Randfold path has been added to \$PATH variable\n"; 
#        `echo 'export PATH=\$PATH:$dir/essentials/randfold-2.0' >> ~/.bashrc`;
#    }

#    $in = `$grep "$dir/essentials/randfold-2.0:*" ~/$shellconf`;
#    if(not $in){ 
#        print STDERR "Randfold path has been added to \$PATH variable\n"; 
#        `echo 'export PATH=\$PATH:$dir/essentials/randfold-2.0' >> ~/$shellconf`;
#}


#    $in2 = `$grep "$dir/essentials/randfold-2.0:*" ~/.cshrc`;
#    if($in2){
	#`echo 'setenv PATH \$PATH:$dir/essentials/randfold-2.0' >> ~/.cshrc`;
	#    }
	chdir "$install_bin_dir";
	if(not -f "randfold"){
		system("ln -s $dir/essentials/randfold-2.0/randfold .");
	}
	chdir "$dir/essentials/"
}


##check for zlib perl
my $zlib=`perl -e 'use Compress::Zlib;' 2>&1`;

if(not $zlib){
	print STDERR "Compress::Zlib\t\t\t\t\t already installed, nothing to do ...\n";
	$progs{zlib}=1;
}else{
	die "please install Compress::Zlib by using CPAN before you proceed\n";

}

#my $pdfapi=`perl -e 'use PDF::API2;' 2>&1`;


$ret = checkBIN("perl -e \'use Font::TTF; print \"installed\";\'","installed");

if($ret == 0){
	print STDERR "Font::TTf                                        already installed, nothing to do ...\n";
	$progs{ttf}=1;
}else{
	my $version='';
	if( -f "CHECKSUMS"){ unlink "CHECKSUMS";}
	`$dtool http://www.cpan.org/authors/id/M/MH/MHOSKEN/CHECKSUMS $dopt`;
	open IN,"CHECKSUMS" or die "File checksums not found\n";
	while(<IN>){
		if(/((Font-TTF-\d.+).tar.gz)/){
			$dfile=$1;
			$version=$2;
		}
	}
	close IN;

	if(not -f $dfile){
		print STDERR "Downloading Font::TTF now\n\n";
		`$dtool http://www.cpan.org/authors/id/M/MH/MHOSKEN/$dfile $dopt`;
	}

	if(not -f $dfile){
		die "Download of Font::TTF failed\n\n";
	}

	print STDERR "Installing Font-TTF now\n\n";
	`tar xzf $dfile`; 
	chdir("$version");

	`perl Makefile.PL INSTALL_BASE=$ENV{'HOME'}/perl5 LIB=$dir/lib/perl5`;
	`make 1>>../install.log 2>>../install_error.log`;

	`mv Makefile Makefile.orig`;

	open IN,"Makefile.orig" or die "No Makefile found\n";
	open OUT,">Makefile" or die "No Makefile found\n";
	while(my $cl= <IN>){
		if($cl =~ /^INSTALL_BASE\s=/){
			print OUT "INSTALL_BASE = $dir\n";
		}else{
			print OUT $cl;
		}
	}
	close IN;

	`make install 1>>../install.log 2>>..install_error.log`;

	$ret = checkBIN("perl -e \'use Font::TTF; print \"installed\";\'","installed");

	if($ret == 0){
		print STDERR "Font::TTF                                        installation successful \n";
		$progs{ttf}=1;
	}	

	chdir("..");
}



$ret = checkBIN("perl -e \'use PDF::API2; print \"installed\";\'","installed");

if($ret == 0){
	print STDERR "PDF::API2                                        already installed, nothing to do ...\n";
	$progs{pdf}=1;
}else{
	my $version='';
	if( -f "CHECKSUMS"){ unlink "CHECKSUMS";}
	`$dtool http://www.cpan.org/authors/id/S/SS/SSIMMS/CHECKSUMS $dopt`;
	open IN,"CHECKSUMS" or die "File checksums not found\n";
	while(<IN>){
		if(/((PDF-API2.+).tar.gz)/){
			$dfile=$1;
			$version=$2;
		}
	}
	close IN;

	if(not -f $dfile){
		print STDERR "Downloading PDF-API2 now\n\n";
		`$dtool http://ftp-stud.hs-esslingen.de/pub/Mirrors/CPAN/authors/id/S/SS/SSIMMS/$dfile $dopt`;
	}

	if(not -f $dfile){
		die "Download of PDF-API2 failed\n\n";
	}

	print STDERR "Installing PDF-API2 now\n\n";
	`tar xzf $dfile`; 
	chdir("$version");

	`perl Makefile.PL INSTALL_BASE=$ENV{'HOME'}/perl5 LIB=$dir/lib/perl5`;
	`make 1>>../install.log 2>>..install_error.log`;
	`mv Makefile Makefile.orig`;

	open IN,"Makefile.orig" or die "No Makefile found\n";
	open OUT,">Makefile" or die "No Makefile found\n";
	while(my $cl= <IN>){
		if($cl =~ /^INSTALL_BASE\s=/){
			print OUT "INSTALL_BASE = $dir\n";
		}else{
			print OUT $cl;
		}
	}
	close IN;

	`make install 1>>../install.log 2>>..install_error.log`;

	$ret = checkBIN("perl -e \'use PDF::API2; print \"installed\";\'","installed");

	if($ret == 0){
		print STDERR "PDF::API2                                        installation successful\n";
		$progs{pdf}=1;
	}

	chdir("..");
}

my $sum=0;
for my $k (keys %progs){
	$sum+=$progs{$k};
	print STDERR "\n\n$k was/is not installed properly\n\n" if(not $progs{$k} and $binnew==0);
}

if($sum == 6){
	print STDERR "\n\nInstallation successful\n\n\n\n\n\n";
	print STDERR "To check if everything works fine you can now change to the 
	tutorial_dir and type 'bash run_tut.sh' to make a test run\n\n";

	chdir $dir;
	open EF,">install_successful" or die "Could not create file install_successful\n
	In case that all tools are running properly then please create this empty file manually 
	in your mirdeep2 installation folder. Otherwise the other tools will not run.
	";
	close EF;
}else{
	print STDERR "\n\nPlease run the install.pl script again to check if 
	everything is properly installed.

	";
}


exit;


sub rem_mirdeep{
	my ($file)=@_;
	if(-f "$ENV{'HOME'}/$file"){
		`cp ~/$file ~/${file}_$time`;
		print STDERR "~/$file backup is ~/${file}_$time\n";
		`mv ~/$file ~/$file.bak`;
		open OUT,">$ENV{'HOME'}/$file" or die "Cannot create file $ENV{'HOME'}/$file\n";
		open IN,"$ENV{'HOME'}/$file.bak";
		my @line;
		my $tmp;
		while(<IN>){
			$tmp="";
			if(/mirdeep/){
				@line = split(/:/);
				foreach(@line){
					if(/mirdeep/){}else{
						$tmp.="$_";
					}
				}
				if($file !~ /.cshrc/){

					if($tmp !~ /PATH=\$PATH$/ and $tmp !~ /PERL5LIB=\$PERL5LIB\s*$/){
						print OUT "$tmp\n";
					}
				}else{
					if($tmp !~ /PATH \$PATH$/ and $tmp !~ /PERL5LIB \$PERL5LIB\s*$/){
						print OUT "$tmp\n";
					}
				}
			}else{
				print OUT;
			}
		}
	}else{
		print "file $ENV{'HOME'}/$file.bak not found\n";
	}
}


sub checkBIN{
	my ($a,$b) = @_;    
	my $e = system("$a> tmp 2>tmp2");

	open IN,"<tmp";
	my $found =1;
	while(<IN>){
		if(/$b/i){
			$found =0;

		}
	}
	close IN;
	if($found){
		open IN,"<tmp2";
		while(<IN>){
			if(/$b/i){
				$found =0;

			}
		}
	}
	close IN;
	return $found;
}


sub check{
	my ($url,$file) = @_;
	my $out='';
	if($dtool =~ /wget/){
		$out=`wget --spider -v $url 2>&1`;
	}elsif($dtool =~ /curl/){
		$out=`curl --head $url |head -n1`;
		if($out =~ /NOT/i){
			$out='broken';
		}

		if($url =~ /https/){
			`curl $url > $file`;
			if(-s "$file" < 1000){
				$out='broken';
				`rm -f $file`;
			}else{
				$out=0;
			}
		}

	}else{
		die "No download tool found\nplease install wget or curl\n";
	}

	if($out =~ /broken/i){
		return 0;
	}else{
		return 1;
	}
}


sub buildgood{
	if(-f $_[0]){
		print STDERR "Building of $_[0] successful\n";
		$progs{$_[1]}=1 if($_[1]);
	}else{
		die "Building of $_[0] not successful\nPlease have a look at the install.log and install_error.log in 
		the essentials directory
		";
	}
}


