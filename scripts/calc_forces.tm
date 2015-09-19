#! /usr/bin/perl
# $Id: calc_forces.tm,v 1.6 2005/06/07 10:12:20 hal Exp $
#
# calc_forces.tm: compute ab initio forces via TURBOMOLE (dscf and grad)
#


use strict;
use subs qw(check_append);


# get script's name
my $progname = $0;
$progname =~ s{^.*/}{};

my $qmmm_coords = $ARGV[0];	# file with coords from QM zone
my $qmmm_forces = $ARGV[1];	# file into which forces will be written
my $qmmm_pntchs = $ARGV[2];	# file containing point charges
my $changed = $ARGV[3];		# flag indicating changes in the QM zone

my $define_exe = "define";	# setup and initial guess
my $dscf_exe = "dscf";		# SCF calculation
my $grad_exe = "grad";		# gradient calculation
my $nprocs = $ENV{'QION_NPROCS'}; # number of processors
my $start_mo = $ENV{'QION_STMO'}; # special method to create start MO's

my $tm_define_in = "tm_define.in"; # must exist before you start!
my $tm_control_file = "control"; # TURBOMOLEs control file
my $tm_coord_file = "coord";	# TURBOMOLEs coordinate file
my $tm_gradient_file = "gradient"; # TURBOMOLEs first derivative file
my $tm_gradient_pntch_file = "gradient_pntch"; # TURBOMOLEs point charge gradients
my $tm_energy_file = "energy"; # TURBOMOLEs energy file
my $tm_append_file = "control.app"; # additional parameters for control file

my $tm_out_file = "tm.out";	# output file for define, dscf and grad
my $tm_out_file_prev = $tm_out_file . ".prev";
my $tm_out = "";		# output of define, dscf, and grad
my $converged_string = 'convergence criteria satisfied after\s*\d+' .
  '\s+iterations';		# search string

my $last_line;
my $energy;

use constant ANG2BOHR => 1.8897259885789232;



# QION_NPROCS may be set to the desired number of processors
$nprocs = 2 unless defined ($nprocs) && $nprocs =~ /^\d$/ && $nprocs > 0;


# transform qmmm coordinates into TURBOMOLE format
open (HOT_COORDS, "< $qmmm_coords") or
  die "$progname: can't open $qmmm_coords";
open (TM_COORDS, "> $tm_coord_file") or
  die "$progname: can't open $tm_coord_file";

print TM_COORDS "\$coord\n";

while (<HOT_COORDS>) {
  my @coords = split;
  printf TM_COORDS ("%20.14f  %20.14f  %20.14f      %s\n",
		    $coords[1] * ANG2BOHR, $coords[2] * ANG2BOHR,
		    $coords[3] * ANG2BOHR, lc ($coords[0]));
}

print TM_COORDS "\$end\n";

close (TM_COORDS);
close (HOT_COORDS);


# check if we have to run define due to changes in the QM zone
if ($changed) {
  unlink "$tm_control_file";

  if ($start_mo eq "dft") {
    my $dft_on = "tm_define_dft.on";
    my $dft_off = "tm_define_dft.off";

    die "$progname: can't find $dft_on" unless -f $dft_on;
    die "$progname: can't find $dft_off" unless -f $dft_off;

    # switch DFT calculation on
    $tm_out = `$define_exe < $dft_on 2>&1 | tee $tm_out_file`;
    exit 1 if ($tm_out !~ /define ended normally/);

    check_append;

    $tm_out = `$dscf_exe $nprocs 2>&1 | tee -a $tm_out_file`;
    exit 1 if ($tm_out !~ /$converged_string/);

    # switch DFT calculation off
    $tm_out = `$define_exe < $dft_off 2>&1 | tee $tm_out_file`;
    exit 1 if ($tm_out !~ /define ended normally/);
  } else {
    die "$progname: can't find $tm_define_in" unless -f $tm_define_in;

    $tm_out = `$define_exe < $tm_define_in 2>&1 | tee $tm_out_file`;
    exit 1 if ($tm_out !~ /define ended normally/);

    check_append;
  }
}

# if point charge file exists append point charges to control file
if (-f $qmmm_pntchs) {
  my $control_contents = "";

  open (CONTROL, "< $tm_control_file") or
    die "$progname: can't open $tm_control_file";

  while (<CONTROL>) {
    last if /^\s*\$point_charges/;
    last if /^\s*\$end/;
    next if /^\s*point charges on/;

    $control_contents .= $_;
    $control_contents .= "   point charges on\n" if /^\s*\$drvopt/;
  }

  close (CONTROL);

  open (CONTROL, "> $tm_control_file") or
    die "$progname: can't open $tm_control_file";
  open (PNT_CHARGE, "< $qmmm_pntchs") or
    die "$progname: can't open $qmmm_pntchs";

  print CONTROL $control_contents;
  print CONTROL "\$point_charges\n";

  while (<PNT_CHARGE>) {
    my @coords = split;
    printf CONTROL ("%14.8f %14.8f %14.8f %8.5f\n",
		    $coords[0] * ANG2BOHR, $coords[1] * ANG2BOHR,
		    $coords[2] * ANG2BOHR, $coords[3]);
  }

  close (PNT_CHARGE);

  print CONTROL "\$point_charge_gradients file=$tm_gradient_pntch_file\n";
  print CONTROL "\$end\n";

  close (CONTROL);
}

# start SCF calculation
$tm_out = `$dscf_exe $nprocs 2>&1 | tee -a $tm_out_file`;

exit 1 if ($tm_out !~ /$converged_string/);


# start gradient calculation
unlink $tm_gradient_file;
unlink $tm_gradient_pntch_file;

$tm_out = `$grad_exe $nprocs 2>&1 | tee -a $tm_out_file`;

exit 1 if ($tm_out !~ /grad ended normally/);


# get energy
open (ENERGY, "< $tm_energy_file") or
  die "$progname: can't open $tm_energy_file\n";
open (FORCES, "> $qmmm_forces") or
  die "$progname: can't open $qmmm_forces";

while (<ENERGY>) {
  next if /\s*\$/;            # ignore control words

  $last_line = $_;
}

close (ENERGY);

# energy in hartree
$energy = (split ' ', $last_line)[1];
print FORCES "$energy\n";

# transform gradient
open (GRADIENT, "< $tm_gradient_file") or
  die "$progname: can't open $tm_gradient_file";

while (<GRADIENT>) {
  next if /^\s*\$/;
  next if /^\s*cycle =/;

  # change exponential D into exponential E
  y/D/E/;

  my @gradient = split;

  if ($#gradient == 2) {
    # write negated gradient
    printf FORCES ("%20.14G  %20.14G  %20.14G\n",
		   -$gradient[0], -$gradient[1], -$gradient[2]);
  }
}

# read point charge gradients if necessary
if (-f $qmmm_pntchs) {
  open (PNTCH, "< $tm_gradient_pntch_file") or
    die "$progname: can't open $tm_gradient_pntch_file";

  while (<PNTCH>) {
    next if (/^\s*\$/);

    y/D/E/;

    my @gradient = split;

    printf FORCES ("%20.14G  %20.14G  %20.14G\n",
		   -$gradient[0], -$gradient[1], -$gradient[2]);
  }

  close (PNTCH);
}

close (FORCES);
close (GRADIENT);

`cp $tm_out_file $tm_out_file_prev`;

# delete temporary files and exit gracefully
unlink $tm_out_file;
unlink $qmmm_pntchs;
unlink $qmmm_coords;

exit 0;


# if append file exists append contents to control file
sub check_append {
  if (-f $tm_append_file) {
    my $contents = "";

    open (CONTROL, "< $tm_control_file") or
      die "$progname: can't open $tm_control_file";

    while (<CONTROL>) {
      last if /^\s*\$end/;

      $contents .= $_;
    }

    close (CONTROL);

    open (CONTROL, "> $tm_control_file") or
      die "$progname: can't open $tm_control_file";
    open (APPEND, "< $tm_append_file") or
      die "$progname: can't open $tm_append_file";

    while (<APPEND>) {
      $contents .= $_;
    }

    close (APPEND);

    print CONTROL "$contents\$end\n";

    close (CONTROL);
  }
}
