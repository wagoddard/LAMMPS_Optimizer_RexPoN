package Packages::BOX;

use FindBin qw($Bin);
use lib "$FindBin::Bin";
require Exporter;
use strict;
use Math::Trig qw(acos);
use Packages::General qw(GetMinDist);
use Cwd;
use constant PI => atan2(1,1) * 4;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(CreateGrid GetNeighbours GetBox PrintBox GetRadii CenterAtoms GetSurface GetNlist InitBox Cart2Frac
			 GetGridDims ConvertBox MakeBox  MoveAtomsToOrigin PlaceAtomsOnGrid GetCellFromHeader box2H Map2UnitCell);
$VERSION = "1.00";

sub getBoxDims {
	my ($BBOX, $grid_len, $printGrid) = @_;
	my ($counter, $bLen);

	if ($printGrid) {
		PrintBox("Total bounding box for atom centers:", $BBOX);
	}

	for $counter (keys %{ $BBOX }) {
		$bLen = $BBOX->{$counter}{"hi"} - $BBOX->{$counter}{"lo"};
		if ($bLen =~ /\d+\.\d+/ || $bLen % $grid_len > 0) {
			$BBOX->{$counter}{"hi"} = ((int($bLen/$grid_len) + 1) * $grid_len) + 
										$BBOX->{$counter}{"lo"};
		}
	}

}

sub CreateGrid {
	my ($ATOMS, $radii, $BBOX, $grid_len, $printGrid) = @_;
	my ($counter, $CURR_BOX, $Index, %GRID, $AtomC, $total_cells, $bLen);
	my ($currCell, $burried, $exclude, $i, $j, $k);

	if (! defined($printGrid)) {
		$printGrid = 0;
	}

	$AtomC = $total_cells = $exclude = $burried = 0;
	if ($radii) {
		getBoxDims($BBOX, $grid_len, $printGrid);
	}
	
	for $i ("X","Y","Z") {
		$Index->{$i} = int($BBOX->{$i}{len}/$grid_len) + 1;
	}
	
	for $i (0 .. $Index->{X}) {
		for $j (0 .. $Index->{Y}) {
			for $k (0 .. $Index->{Z}) {
				$total_cells++;
				$GRID{$i}{$j}{$k} = (
									{
										"XINDEX" => $i,
										"YINDEX" => $j,
										"ZINDEX" => $k,
										"X"	  => (
												{
													"lo" => ($BBOX->{X}{lo} + ($i * $grid_len)),
													"hi" => ($BBOX->{X}{lo} + (($i+1) * $grid_len)),
												}
												),
										"Y"	  => (
												{
													"lo" => ($BBOX->{Y}{lo} + ($j * $grid_len)),
													"hi" => ($BBOX->{Y}{lo} + (($j+1) * $grid_len)),
												}
												),
										"Z"	  => (
												{
													"lo" => ($BBOX->{Z}{lo} + ($k * $grid_len)),
													"hi" => ($BBOX->{Z}{lo} + (($k+1) * $grid_len)),
												}
												)
									}
									);
			}
		}
	}

	($GRID{X}{tot}, $GRID{Y}{tot}, $GRID{Z}{tot}) = ($Index->{X}, $Index->{Y}, $Index->{Z});
	($AtomC, $exclude) = PlaceAtomsOnGrid($ATOMS, \%GRID, $BBOX, $grid_len) if (defined($ATOMS));

	if ($printGrid) {
		PrintBox("Total vdw box size:", $BBOX);
		print "Total Atoms: $AtomC\n";
		print "Total Cells: $total_cells\n";
		print "Cells with Atoms: $exclude (" . sprintf("%.2f", (100 * ($exclude/$total_cells))) . ")%\n";
	}
	return (\%GRID, $BBOX, $AtomC);
}

sub setAllBurriedCells {
	my ($GRID) = $_[0];
	my ($x, $y, $z, $cell, $burried, $CLIST, $nCells, $newBurried, $j, $newCell);
	my ($xIndex, $yIndex, $zIndex);

	$burried = 0;
	$newBurried = 1;
	while ($newBurried) {
		$newBurried = 0;
		for $x (keys %{ $GRID }) {
			for $y (keys %{ $GRID->{$x} }) {
				CELLLOOP: for $z (keys %{ $GRID->{$x}{$y} }) {
					$cell = \%{ $GRID->{$x}{$y}{$z} };
					if ($cell->{Burried}) {
						$cell->{Surface} = 0;
						$burried++;
						next CELLLOOP;
					}
					$CLIST = GetNeighbours($GRID, $cell);
					for $nCells (@{ $CLIST }) {
						$xIndex = $nCells->{XINDEX};
						$yIndex = $nCells->{YINDEX};
						$zIndex = $nCells->{ZINDEX};
						$newCell = \%{ $GRID->{$nCells->{XINDEX}}{$nCells->{YINDEX}}{$nCells->{ZINDEX}} };
						next CELLLOOP if (! exists($newCell->{VOL}) || $newCell->{VOL} == 0 || ! $newCell->{Burried});
					}
					$cell->{Burried} = 1; #if all of my neighbours are burried then i am as well
					$cell->{Surface} = 0;
					$burried++;
					$newBurried = 1;
				}
			}
		}
	}
	return $burried;
}
				
sub PlaceAtomsOnGrid {
	my ($ATOMS, $GRID, $BOX, $grid_len) = @_;
	my ($counter, $AtomC, $Index, $dim, @deleted_keys, $i);  
	my ($atom, $cell, $maxVol, $exclude, $atomVol, $burried);
	my ($pi) = atan2(1,1) *4;
	
	$maxVol = $grid_len**3;

	$AtomC = $counter = $exclude = $burried = 0;
	@deleted_keys = keys %{ $ATOMS };
	while ($#deleted_keys > -1) {
		$counter = pop @deleted_keys;
		for $dim ("X", "Y", "Z") {
			$i = int(($ATOMS->{$counter}{$dim . "COORD"} - $BOX->{$dim}{lo})/$grid_len);
			$Index->{$dim} = $i;
		}
		$AtomC++;
		$atom = \%{ $ATOMS->{$counter} };
		next if (! exists($GRID->{$Index->{X}}) or
				 ! exists($GRID->{$Index->{X}}{$Index->{Y}}) or
				 ! exists($GRID->{$Index->{X}}{$Index->{Y}}{$Index->{Z}}));
		$cell = \%{ $GRID->{$Index->{X}}{$Index->{Y}}{$Index->{Z}} };
		$atom->{CELL}{XINDEX} = $Index->{X};
		$atom->{CELL}{YINDEX} = $Index->{Y};
		$atom->{CELL}{ZINDEX} = $Index->{Z};
		
		$atomVol  = 0;
		$cell->{VOL} = 0 if (! exists($cell->{VOL}));
		if ($atom->{IS_SOLVENT}) {
			push @{ $cell->{WATERS} }, $ATOMS->{$counter};
		} elsif ($atom->{NUMBONDS} == 0) {
			push @{ $cell->{IONS} }, $ATOMS->{$counter};
		} elsif ($atom->{RESNAME} !~ /(WAT|Na|Mg|Cl)/i) {
			push @{ $cell->{ATOMS} }, $ATOMS->{$counter};
			$atomVol  = (4 * $atom->{"RADII"}**3 * $pi)/3;
			$exclude++ if ($#{ $cell->{ATOMS} } == 0);
		} elsif ($atom->{RESNAME} eq "WAT" or $atom->{RESNAME} eq "HOH") {
			push @{ $cell->{WATERS} }, $ATOMS->{$counter};
		} else {
			push @{ $cell->{IONS} }, $ATOMS->{$counter};
		}
		$cell->{VOL} += $atomVol;
		&setExcludedVol($cell) if ($atomVol > 0);
	}
	return ($AtomC, $exclude);
}

sub setExcludedVol {
	my ($cell) = $_[0];
	my ($atom, $radii, $dim);

	for $atom (@{ $cell->{ATOMS} }, @{ $cell->{IONS} }) {
		$radii = $atom->{RADII};
		for $dim ("X", "Y", "Z") {
			$cell->{EXCLUDE}{$dim}{hi} = ($atom->{$dim . "COORD"} + $radii)
				if (! exists($cell->{EXCLUDE}{$dim}{hi}) or 
					($atom->{$dim . "COORD"} + $radii) > $cell->{EXCLUDE}{$dim}{hi});
			$cell->{EXCLUDE}{$dim}{"lo"} = ($atom->{$dim . "COORD"} - $radii)
				if (! exists($cell->{EXCLUDE}{$dim}{lo}) or 
					($atom->{$dim . "COORD"} - $radii) < $cell->{EXCLUDE}{$dim}{lo});
		}
	}
}

sub GetNeighbours {
	my ($GRID, $CELL, $layers, @FROZEN) = @_;
	my (@CLIST, $i, %MAPING, $j, $count, $CONSTANT);
	my ($xI, $yI, $zI, $k, $tot, $VALS);
	
	$CONSTANT = () if (! @FROZEN);
	for $i (@FROZEN) {
		$CONSTANT->{$i} = 1;
	}
	$layers = 1 if (! defined($layers));
	for $i ("XINDEX", "YINDEX", "ZINDEX") {
		return () if (! exists($CELL->{$i}));
	}

	($xI, $yI, $zI) = ($CELL->{XINDEX}, $CELL->{YINDEX}, $CELL->{ZINDEX});

	for $i ("X", "Y", "Z") {
		$VALS->{$i}{TOT} = eval('$GRID->{' . $i . '}{tot}');
		$VALS->{$i}{INDEX} =  eval('$' . lc($i) . "I");
		if (exists($CONSTANT->{$i})) {
			push @{ $MAPING{$i} }, $VALS->{$i}{INDEX};
			next;
		}
		for $j ((-1 * $layers) .. $layers) {
			$k = $j + $VALS->{$i}{INDEX};
			$tot = $VALS->{$i}{INDEX};
			if ($k < 1) {
				$k += $tot;
			} elsif ($k > $tot) {
				$k -= $tot;
			}
			push @{ $MAPING{$i} }, $k;
		}
	}

	for $i (@{ $MAPING{X} }) {
		for $j (@{ $MAPING{Y} }) {
			for $k (@{ $MAPING{Z} }) {
				push @CLIST, $GRID->{$i}{$j}{$k};
			}
		}
	}
	
	return \@CLIST;
}

sub IsInBox(@) {
	my ($Atom, $Box) = @_;
	my ($index, $returnval);

	$returnval = 0;
	if ($Atom->{"XCOORD"} >= $Box->{"x1"} and $Atom->{"XCOORD"} <= $Box->{"x2"}) {
		if ($Atom->{"YCOORD"} >= $Box->{"y1"} and $Atom->{"YCOORD"} <= $Box->{"y2"}) {
			if ($Atom->{"ZCOORD"} >= $Box->{"z1"} and $Atom->{"ZCOORD"} <= $Box->{"z2"}) {
				$returnval = 1;
			}
		}
	}

	return $returnval;

}

sub Cart2Frac {
	my ($atoms, $box) = @_;
	my ($i, $ori, $F);

	$F = $box->{F};
	for $i (keys %{ $atoms }) {
		$ori->[0] = $atoms->{$i}{XCOORD};
		$ori->[1] = $atoms->{$i}{YCOORD};
		$ori->[2] = $atoms->{$i}{ZCOORD};
		$atoms->{$i}{FA} = $ori->[0]*$F->[0][0]+$ori->[1]*$F->[0][1]+$ori->[2]*$F->[0][2];
		$atoms->{$i}{FB} = $ori->[0]*$F->[1][0]+$ori->[1]*$F->[1][1]+$ori->[2]*$F->[1][2];
		$atoms->{$i}{FC} = $ori->[0]*$F->[2][0]+$ori->[1]*$F->[2][1]+$ori->[2]*$F->[2][2];
	}
}

sub InitBox {
	my ($box, $atoms) = @_;

	&box2H($box);
	&box2F($box);
	&Cart2Frac($atoms, $box);
}

sub GetBox {
	my ($ATOMS, $PARMS, $HEADERS) = @_;
	my (%BOX, $isValid, $i, $j, $coord); 
	my ($cVal, @tmp, $radii, @vals);
	
	$isValid = 0;
	@tmp = (-9999,99999,90,-99999,99999,90,-99999,99999,90);
	if (defined($HEADERS)) {
		for $i (@{ $HEADERS }) {
			if ($i =~ /^CRYSTX\s+(.+)/) {
				@vals = split /\s+/, $1;
				for $j (0 .. 2) {
					$tmp[($j * 3)] = $vals[$j];
					$tmp[($j * 3) + 1] = 0,
					$tmp[($j * 3) + 2] = $vals[3+$j];
				}
				$isValid = 1;
			}
		}
	}
	
	if (! $isValid) {
		for $i (keys %{ $ATOMS }) {
			$radii = GetRadii($ATOMS->{$i}, $PARMS);
			if (! defined($radii)) {
				$radii = 0;
			}
			
			$j = 0;
			for $coord ("X", "Y", "Z") {
				$cVal = $ATOMS->{$i}{$coord . "COORD"};
				$tmp[$j] = ($cVal + $radii)
					if (($cVal + $radii) > $tmp[$j]);
				$tmp[$j + 1] = ($cVal - $radii)
					if (($cVal - $radii) < $tmp[$j + 1]);
				$j += 3;
			}
		}
	}

	%BOX = (
			"X" => {
				"hi"	=> $tmp[0],
				"lo"	=> $tmp[1],
				"len"   => $tmp[0] - $tmp[1],
				"angle" => $tmp[2],
			},
			"Y" => {
				"hi"	=> $tmp[3],
				"lo"	=> $tmp[4],
				"len"   => $tmp[3] - $tmp[4],
				"angle" => $tmp[5],
			},
			"Z" => {
				"hi"	=> $tmp[6],
				"lo"	=> $tmp[7],
				"len"   => $tmp[6] - $tmp[7],
				"angle" => $tmp[8],
			},
		);
	&InitBox(\%BOX, $ATOMS);	
	return \%BOX;
}

sub GetRadii {
	my ($atom, $PAR) = @_;
	my ($returnVal, $atmName);

	$atmName = $atom->{"FFTYPE"};
	if (! $PAR or ! exists($PAR->{"VDW"}{$atmName}) or ! exists($PAR->{"VDW"}{$atmName}{$atmName})) {
		$atmName = $atom->{"ATMNAME"};
	} else {
		$returnVal = $PAR->{"VDW"}{$atmName}{$atmName}{1}{"VALS"}[1]/2;
	}

	return $returnVal;
}

sub PrintBox {
	my ($intext, $inBox) = @_;
	my ($counter, $out_string);

	$out_string = sprintf("%-40s", $intext);
	for $counter ("X", "Y", "Z") {
		$out_string .= sprintf("%8.3f ", ($inBox->{$counter}{"hi"} - $inBox->{$counter}{"lo"}));
	}

	print "$out_string\n";

}

sub CenterAtoms {
	my ($ATOMS, $BOX, $radii) = @_;
	my (%Offset, $dim, $atomC);

	for $dim ("X", "Y", "Z") {
		$Offset{$dim} = $BOX->{$dim}{"lo"};
	}

	for $atomC (keys %{ $ATOMS }) {
		for $dim ("X", "Y", "Z") {
			$ATOMS->{$atomC}{$dim . "COORD"} -= $Offset{$dim};
		}
	}
   
}

sub GetSurface {
	my ($GRID) = $_[0];
	my ($i, $j, $k,%SURFACE, $sCell, $cellVol, $atom);
	my ($NEIGH, %SCHECK, $currCell, $count, $l, $dim);  
	my ($bCell, $charge, $totCharge, $CELLS, $ions);

	$cellVol = 1;
	$bCell = $ions = $sCell = 0;
	for $i ("X", "Y", "Z") {
		$cellVol *= ($GRID->{1}{1}{1}{$i}{hi}-$GRID->{1}{1}{1}{$i}{lo});
	}
	$sCell = $charge = $totCharge = $ions = 0;
	for $i (1 .. $GRID->{X}{tot}) {
		for $j (keys %{ $GRID->{$i} }) {
			for $k (keys %{ $GRID->{$i}{$j} }) {
				if (exists($GRID->{$i}{$j}{$k}{ATOMS})) {
					$GRID->{$i}{$j}{$k}{Surface} = 1;
					$SCHECK{$i}{$j}{$k} = $GRID->{$i}{$j}{$k};
				}
			}
		}
	}

	for $i (keys %SCHECK) {
		for $j (keys %{ $SCHECK{$i} }) {
			for $k (keys %{ $SCHECK{$i}{$j} }) {
				$currCell = $GRID->{$i}{$j}{$k};
				$CELLS = GetNeighbours($GRID, $currCell);
				$count = 0;
				for $l (@{ $CELLS }) {
					$count++ if (exists($l->{Surface}));
				}
				$count--;
				if ($count == 26) {
					$bCell++;
					delete($currCell->{Surface});
					$currCell->{Burried} = 1;
					next;
				}
				$sCell++;
				$charge = 0;
				for $dim ("X", "Y", "Z") {
					$SURFACE{$sCell}{$dim} = $currCell->{$dim . "INDEX"};
				}
				for $atom (@{ $currCell->{ATOMS} }, @{ $currCell->{IONS} }) {
					$charge += $atom->{CHARGE};
				}
				$SURFACE{$sCell}{CHARGE} = $charge;
				$SURFACE{$sCell}{BURRIED} = 0;
				$totCharge += $charge;
				$ions++ if (exists($currCell->{IONS}));
			}
		}
	}

	print "Burried Cells: $bCell\n";
	print "Cells defining molecular surface: $sCell\n";
	print "Surface Cells containing ions: $ions\n";
	printf "Total charge of surface atoms: %.3f\n", $totCharge;

	return \%SURFACE;

}

sub getTriclinicParms {
	my ($box, $cell) = @_;
	my ($hyp, $angle);
	my ($PI) = atan2(1,1) * 4;

	$angle = $box->{$cell->{3}}{angle};
	$angle = $angle % 180;
	$angle = 180 - $angle if ($angle > 90);
	$angle *= ($PI/180);
	$hyp = $box->{$cell->{2}}{len};

	return ($angle, $hyp);
}

sub GetGridDims {
	my ($molBox, $dimOpt, $min, $max) = @_;
	my ($i, %GRID, $j, $angle, $c1, $c2, $hyp, @TRICLINIC);

	@TRICLINIC = (
					(
					{
						1 => "X",
						2 => "Y",
						3 => "Z",
					},
					),
					(
					{
						1 => "X",
						2 => "Z",
						3 => "Y",
					},
					),
					(
					{
						1 => "Y",
						2 => "Z",
						3 => "X",
					},
					)
					);

#corrections for triclinic cell
	for $i (@TRICLINIC) {
		($angle, $hyp) = getTriclinicParms($molBox,$i);
		next if (cos($angle) < 0.01);
		$c1 = $hyp * cos($angle)/2;
		$c2 = $hyp * (1 - sin($angle)); 
		$molBox->{$i->{2}}{lo} += $c1;
		$molBox->{$i->{2}}{hi} -= $c1;
		$molBox->{$i->{2}}{len} -= ($c1 * 2);

		$molBox->{$i->{1}}{hi} -= $c2;
		$molBox->{$i->{1}}{lo} += $c2;
		$molBox->{$i->{1}}{len} -= ($c2 * 2);
	}

	for $i (keys %{ $molBox }) {
		for $j ("hi", "lo") {
			$GRID{$i}{$j} = $molBox->{$i}{$j};
		}
	}

	for $i ("x","y","z") {
		$j = uc($i);
		if (lc($dimOpt) =~ /$i/) {
			if ($min =~ /\+(\d+\.*\d*)/) {
				$GRID{$j}{lo} -= $1;
			}elsif ($min < $molBox->{$j}{lo}) {
				$GRID{$j}{lo} = $min;
			}

			if ($max =~ /\+(\d+\.*\d*)/) {
				$GRID{$j}{hi} += $1;
			}elsif ($max > $molBox->{$j}{hi}) {
				$GRID{$j}{hi} = $max;
			}
		}
	}

	return \%GRID;
}

sub MakeBox { 
	my ($box) = $_[0];
	my (%BOX, $i, @dim);

	@dim = ("X", "Y", "Z");
	for $i (0 .. 2) {
		$BOX{$dim[$i]}{lo} = 0;
		$BOX{$dim[$i]}{hi} = $box->{$i + 2}{DATA};
		$BOX{$dim[$i] . "COORD"}{lo} = 0;
		$BOX{$dim[$i] . "COORD"}{hi} = $box->{$i + 2}{DATA};
		$BOX{$dim[$i] . "COORD"}{len} = $box->{$i + 2}{DATA};
	}
	$BOX{ALPHA}=$BOX{BETA}=$BOX{GAMMA} = 90;
	return \%BOX;
}

sub ConvertBox {
	my ($lammpsBox) = $_[0];
	my (%amberBox, $i);

	$amberBox{1}{DATA} = 90;
	for $i (0 .. 2) {
		$amberBox{$i + 2}{DATA} = $lammpsBox->[$i]{hi} - $lammpsBox->[$i]{lo};
	}

	return \%amberBox;
}

sub MoveAtomsToOrigin {
	my ($atoms) = $_[0];
	my ($MIN, $i, $j);

	for my $i (keys %{ $atoms }) {
		for $j ("XCOORD", "YCOORD", "ZCOORD") {
			$MIN->{$j} = $atoms->{$i}{$j} if (! exists($MIN->{$j}) or $atoms->{$i}{$j} < $MIN->{$j});
		}
	}

	for my $i (keys %{ $atoms }) {
		for $j ("XCOORD", "YCOORD", "ZCOORD") {
			$atoms->{$i}{$j} -= $MIN->{$j};
		}
	}

}

sub GetCellFromHeader {
	my ($headers) = $_[0];
	my ($box, $i);


	for $i (@{ $headers }) {
		(undef, $box->{X},$box->{Y},$box->{Z},$box->{alpha},$box->{beta},$box->{gamma}) = split /\s+/,$i
			if ($i =~ /CRYSTX/);
	}

	die "ERROR: No box info found in bgf file!\n" if (!defined($box));

	return $box;
}

sub box2F {
	my ($box) = $_[0];
	my ($la,$lb,$lc,$alpha,$beta,$gamma,$vol);

	$la=$box->{X}{len}; $lb=$box->{Y}{len}; $lc=$box->{Z}{len};
	$alpha=$box->{X}{angle}*PI/180; $beta=$box->{Y}{angle}*PI/180; $gamma=$box->{Z}{angle}*PI/180;
	$vol = $la*$lb*$lc*sqrt(1-cos($alpha)**2-cos($beta)**2-cos($gamma)**2+2*cos($alpha)*cos($beta)*cos($gamma));

	$box->{F}[0][0] = 1/$la; 
	$box->{F}[0][1] = -cos($gamma)/$la/sin($gamma); 
	$box->{F}[0][2] = (cos($alpha)*cos($gamma)-cos($beta))/$la/$vol/sin($gamma);
	$box->{F}[1][0] = 0;
	$box->{F}[1][1] = 1/$lb/sin($gamma);
	$box->{F}[1][2] = (cos($beta)*cos($gamma)-cos($alpha))/$lb/$vol/sin($gamma);
	$box->{F}[2][0] = $box->{F}[2][1] = 0;
	$box->{F}[2][2] = $la*$lb*sin($gamma)/$vol;

	$box->{Finv}[0][0] = $la; 
	$box->{Finv}[0][1] = $lb*cos($gamma);
	$box->{Finv}[0][2] = $lc*cos($beta);
	$box->{Finv}[1][0] = 0;
	$box->{Finv}[1][1] = $lb*sin($gamma);
	$box->{Finv}[1][2] = $lc*(cos($alpha)-cos($beta)*cos($gamma))/sin($gamma);
	$box->{Finv}[2][0] = $box->{Finv}[2][1] = 0;
	$box->{Finv}[2][2] = $vol/$la/$lb/sin($gamma);
}


sub box2H {
	my ($box) = $_[0];
	my ($aa,$bb,$cc,$fac);
	my (@a,@b,@c,$la,$lb,$lc,$alpha,$beta,$gamma);

	$la=$box->{X}{len}; $lb=$box->{Y}{len}; $lc=$box->{Z}{len};
	$alpha=$box->{X}{angle}; $beta=$box->{Y}{angle}; $gamma=$box->{Z}{angle};
	
	$fac=acos(-1.0)/180.0;
	$aa=$alpha*$fac; $bb=$beta*$fac; $cc=$gamma*$fac;
	#place c along Z
	$c[0]=0; $c[1]=0; $c[2]=$lc;
	#place b in the yz plane
	$b[0]=0; $b[1]=$lb*sin($aa); $b[2]=$lb*cos($aa);
	#find vector for a
	$a[2]=$la*cos($bb);
	$a[1]=($la*$lb*cos($cc)-$a[2]*$b[2])/$b[1];  #la*(cos(cc)-cos(aa)*cos(bb))/sin(aa);
	$a[0]=sqrt($la*$la-$a[1]*$a[1]-$a[2]*$a[2]);

	$box->{H}[0][0]=$a[0]; $box->{H}[1][0]=$a[1]; $box->{H}[2][0]=$a[2];
	$box->{H}[0][1]=$b[0]; $box->{H}[1][1]=$b[1]; $box->{H}[2][1]=$b[2];
	$box->{H}[0][2]=$c[0]; $box->{H}[1][2]=$c[1]; $box->{H}[2][2]=$c[2];

	$box->{Hinv}[0][0]=1/$box->{H}[0][0];
	$box->{Hinv}[1][1]=1/$box->{H}[1][1];
	$box->{Hinv}[2][2]=1/$box->{H}[2][2];

	$box->{Hinv}[1][0]=-$box->{H}[1][0]/$box->{H}[0][0]/$box->{H}[1][1];
	$box->{Hinv}[2][0]=($box->{H}[1][0]*$box->{H}[2][1]-$box->{H}[1][1]*$box->{H}[2][0])/$box->{H}[0][0]/$box->{H}[1][1]/$box->{H}[2][2];
	$box->{Hinv}[2][1]=-$box->{H}[2][1]/$box->{H}[1][1]/$box->{H}[2][2];
	$box->{Hinv}[0][1]=0;
	$box->{Hinv}[1][2]=0;
	$box->{Hinv}[0][2]=0;

}

sub Map2UnitCell {
	my ($coord, $box) = @_;
	my ($Hinv, $H, $n, $ori);

	$Hinv = $box->{Hinv};
	$H = $box->{H};
	@{ $ori }  = ($coord->{XCOORD}, $coord->{YCOORD}, $coord->{ZCOORD});

	$n->[0]=$Hinv->[0][0]*$ori->[0]+$Hinv->[0][1]*$ori->[1]+$Hinv->[0][2]*$ori->[2];
	$n->[1]=$Hinv->[1][0]*$ori->[0]+$Hinv->[1][1]*$ori->[1]+$Hinv->[1][2]*$ori->[2];
	$n->[2]=$Hinv->[2][0]*$ori->[0]+$Hinv->[2][1]*$ori->[1]+$Hinv->[2][2]*$ori->[2];

	while($n->[0]>1) { $n->[0]-=1; }
	while($n->[0]<0) { $n->[0]+=1; }
	while($n->[1]>1) { $n->[1]-=1; }
	while($n->[1]<0) { $n->[1]+=1; }
	while($n->[2]>1) { $n->[2]-=1; }
	while($n->[2]<0) { $n->[2]+=1; }

	$coord->{SHIFT}{XCOORD}=$H->[0][0]*$n->[0]+$H->[0][1]*$n->[1]+$H->[0][2]*$n->[2]-$ori->[0];
	$coord->{SHIFT}{YCOORD}=$H->[1][0]*$n->[0]+$H->[1][1]*$n->[1]+$H->[1][2]*$n->[2]-$ori->[1];
	$coord->{SHIFT}{ZCOORD}=$H->[2][0]*$n->[0]+$H->[2][1]*$n->[1]+$H->[2][2]*$n->[2]-$ori->[2];

	$coord->{XCOORD} += $coord->{SHIFT}{XCOORD};
	$coord->{YCOORD} += $coord->{SHIFT}{YCOORD};
	$coord->{ZCOORD} += $coord->{SHIFT}{ZCOORD};
}

sub GetNlist {
	my ($atoms, $select, $box, $cut, $nlist) = @_;
	my ($i, $j, $k, $dist, $boffset);
	my ($atom1, $atom2, $rec);

	for $i (keys %{ $select }) {
		$atom1 = $atoms->{$i};
		for $j (keys %{ $select }) {
			next if ($i <= $j);
			$atom2 = $atoms->{$j};
			$rec = ();
			&GetMinDist($atom1, $atom2, $box, \%{ $rec });
			$dist = 0.0;
			for $k ("X","Y","Z") {
				$dist += $rec->{$k}*$rec->{$k};
			}
			$dist = sqrt($dist);
			$rec->{r} = $dist;
			next if ($dist>$cut); 
			$rec->{theta} = acos($rec->{Z}/$rec->{r});
			$rec->{phi} = 0;
			$rec->{phi} = atan2($rec->{Y}/$rec->{X},1) if ($rec->{X} != 0);
			$nlist->{$i}{$j} = $rec;
			$rec->{theta} = acos(-$rec->{Z}/$rec->{r});
			$nlist->{$j}{$i} = $rec;
		}
	}
}

1;
