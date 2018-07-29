package Packages::General;

require Exporter;
use strict;
use Packages::CERIUS2 qw(parseCerius2FF parseReaxFF);
use Packages::MPSIM qw(parseMPSimFF);
use Storable qw(dclone);
use FindBin qw($Bin);
use lib "$FindBin::Bin";
use POSIX "fmod";

use Cwd;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(IsInteger Swap CenterText STDev IsDecimal IsBoolean Trim ParseSanderInfoFile TrjSelections
				GetMPSimEnergy FileTester LoadElements val CombineMols CreateGrid GetBondLength GetAngle 
				execCmd CenterOnMol GetTorsion CrossProduct DotProduct CRadDegrees Rotate PrintProgress 
				CoP GetSelections CoM ParseSelect GetSoluteAtoms Permutate ShowSelectionInfo LoadFFs 
				GetStats GetEquilPoint GetTime Round GetSigDigits getFFtype ReadParmFile printAtomSelectOptions
				GetRotations GetRotationMatrix ShuffleArray VecLen Normalize GetMinDist dPcmplx HasCell
				AddElementField FindElement AddRingField);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub GetSigDigits {
	my ($inStr) = $_[0];
	my ($isDec, $count);

	$inStr += 0;
	$isDec = 0;
	if ($inStr =~ /\.(.+)$/) {
		$inStr = "$1";
		$isDec = 1;
	}
	
	$count = length($inStr);
	while ($count > 0) {
		if (substr($inStr, ($count - 1), 1) eq "0") {
			  $count--;
		} else {
			last;
		}
	}

	$count *= -1 if (! $isDec);
	return $count;
}

sub Round {
	my ($n, $places) = @_;
	my $factor = 10 ** ($places || 0);
	return int(($n * $factor) + ($n < 0 ? -1 : 1) * 0.5) / $factor;
}

sub GetStats {
	my ($data, $start) = ($_[0], $_[1]);
	my (%stats, $i, $numPoints, $sigma, $x1, $x2);

	return () if (! @{ $data });

	$start = 1 if (! $start);
	$start--;
	$numPoints = 0;
	for $i ($start .. $#{ $data }) {
		$x1 += $data->[$i];
		$x2 += ($data->[$i]*$data->[$i]);
		$numPoints++;
	}

	$x1 /= $numPoints;
	$x2 /= $numPoints;
	$sigma = ($x2 - $x1*$x1);
	if ($sigma > 0.00001)  {
		$sigma = sqrt($sigma);
	} else {
		$sigma = 0;
	}
	$stats{STDEV} = $sigma;
	$stats{AVG} = $x1;
	$stats{NUM} = $numPoints;
	$stats{TOT} = $stats{AVG} * $stats{NUM};
	
	return \%stats;
}

sub GetEquilPoint {
	my ($data, $tol, $hashKey) = @_;
	my (@sortedData, $finalAvg, $tot, $i, $window, %residuals, $res, $dev);
	my ($equilPoint, @tmp, $j, $isFound);

	$tot = scalar keys %{ $data };
	@tmp = sort numerically keys %{ $data };
	for $i (0 .. ($tot -1)) {
		if (! defined($hashKey)) {
			$sortedData[$i] = $data->{$tmp[$i]};
		} else {
			$sortedData[$i] = $data->{$tmp[$i]}{$hashKey};
		}
	}

	$finalAvg = 0;
	if ($tot < 10) {
		$finalAvg = $sortedData[$tot];
	} else {
		$window = 5;
		if ($tot > 100) {
			$window += sprintf("%.0f", ($tol * $tot));
		}
		for $i (($tot - $window - 1) .. ($tot - 1)) {
			$finalAvg += $sortedData[$i];
		}
		$finalAvg /= $window;
	}
		
	$dev = $tol * $finalAvg;
	
	for $i (0 .. ($tot -1)) {
		$res = abs($sortedData[$i] - $finalAvg);
		$residuals{$i} = $res if ($res < $dev);
	}

	if (! %residuals) {
		return (0, \@sortedData);
	}
	@tmp = reverse sort numerically keys %residuals;
	$equilPoint = pop @tmp;
	
	for $i (@tmp) {
		if (($equilPoint - $i) > 1) {
			$isFound = 0;
			for $j (1 .. $window) {
				if (($equilPoint - $i - $j) == 1) { # found val with window so continue
					$isFound = 1;
					last;
				}
			}
			last if (! $isFound); #if not found in window, equil point is found
		}
		$equilPoint = $i;
	}

	for $i (0 .. ($equilPoint - 1)) {
		splice @sortedData, 0, 1;
	}
   
	$equilPoint = 0 if (! $isFound);
	return ($equilPoint, \@sortedData);
}

sub getFFtype {
	my ($ffFile) = $_[0];
	my ($ffType) = "CERIUS2";
	my ($testCmd) = "grep -c CERIUS $ffFile";

	open TESTCMD, "$testCmd |" or die "ERROR: Cannot execute $testCmd: $!\n";
	while (<TESTCMD>) {
		chomp;
		if ($_ =~ /(\d+)/) {
			$ffType = "MPSIM" if (! $1);
			last;
		}
	}
	return $ffType;
}
	
sub LoadFFs {
	my ($ffs, $step, $alter) = @_;
	my ($i, $PARMS, $counter, $len, $printStr, $ffType, $ELEMENTS);
																														 
	$PARMS = ();
	$alter = 1 if (! defined($alter) or $alter !~ /^1|0$/);
	$ELEMENTS = &LoadElements;
	$counter = 0;
	for $i (@{ $ffs }) {
		print "Step $step: " if (defined($step));
		if (ref($i) ne "HASH") {
			$ffType = getFFtype($i);
			print "Loading $ffType force field $i...";
			$len = 43 + length($i);
			if ($ffType eq "CERIUS2") {
				$PARMS = parseCerius2FF($i, $alter, $PARMS);
			} elsif ($ffType eq "MPSIM") {
				$PARMS = parseMPSimFF($i, $ELEMENTS, $alter, $PARMS);
			} elsif ($ffType eq "REAX") {
				$PARMS = parseReaxFF($i, $ELEMENTS, $alter, $PARMS);
			}
		} else {
			print "Loading $i->{FFTYPE} force field $i->{FF}...";
			$len = 43 + length($i->{FF});
			if ($i->{FFTYPE} eq "CERIUS2") {
				$PARMS = &parseCerius2FF($i->{FF},$alter, $PARMS);
			} elsif ($i->{FFTYPE} eq "MPSIM") {
				$PARMS = &parseMPSimFF($i->{FF}, $ELEMENTS, $alter, $PARMS);
			} elsif ($i->{FFTYPE} eq "REAX") {
				$PARMS = &parseReaxFF($i->{FF}, $ELEMENTS, $alter, $PARMS);
			}
			if (! defined($ffType)) {
				$ffType = $i->{FFTYPE};
			} elsif ($ffType ne "CERIUS2") {
				$ffType = "MULTIPLE";
			}
																														 
		}
		$counter++;
		print "Done\r";
	}
	die "ERROR: No valid $ffType forcefields found!\n" if (! keys %{ $PARMS });
	$printStr = "Step $step: " if (defined($step));
	$printStr .= "Loading $ffType force field..sucessfully loaded $counter force field";
	$printStr .= "s" if ($counter > 1);
	print "$printStr";
	printf("%" . ($len - length($printStr) + 5) . "s", " ") if ($len > length($printStr));
	print "\n";
																														 
	return $PARMS;
}

sub Trim {
	my ($inString) = $_[0];
																				
	for ($inString) {
		s/^\s+//;
		s/\s+$//;
	}
																				
	return $inString;
}

sub FileTester {
	my ($inFile) = $_[0];

	die "Error accessing regular file $inFile: $!\n"
		if (! -e $inFile or ! -r $inFile or ! -T $inFile);
	
}

																			 
sub IsBoolean {

	my ($inString) = $_[0];
	
	$inString =~ /^[1|0|true|false]/i ? 
		return 1 :
		return 0;

}

sub IsDecimal {

	my ($inString) = $_[0];

	return 0 if ($inString !~ /^\-?\d+\.(\d+)$/);
	return 0 if (! $1);
	return 1;
}

sub IsInteger {
	my ($inString) = $_[0];
	
	$inString =~ /^\-?\d+$/ ? 
		return 1 :
		return 0;
}

sub Swap {
	my ($val1, $val2) = @_;
	return ($val2, $val1);
}

sub CenterText {
	my ($in_text, $field_lenght) = @_;
	my ($text_length, $result);

	$text_length = $field_lenght - length($in_text) ;

	if ($text_length > 0) {
		$result = sprintf("%" . ( int($text_length/2) + ($text_length % 2) ) . "s", " ");
		$result .= sprintf("%" . length($in_text) . "s", "$in_text");
		$result .= sprintf("%" . int($text_length/2) . "s", " ");
	}else {
		$result = $in_text;
	}

	return $result;
}

sub STDev {
	my ($dataString) = $_[0];
	my ($x1, $x2, $i, $count, $sigma);

	$count = $sigma = 0;
	#while ($dataString =~ /(\-?\d+\.?\d*e?\-?\+?\d*)/g) {
	while ($dataString =~ /(\S+)/g) {
		$count++;
		$x1 += $1;
		$x2 += ($1*$1);		
	}

	if ($count > 0) {
		$x1 /= $count;
		$x2 /= $count;
		$sigma = ($x2 - $x1*$x1);
	}

	if ($sigma > 0.00001) {
		$sigma = sqrt($sigma);
	} else {
		$sigma = 0;
	}
	return ($x1, $sigma, ($x1 * $count));
}

sub STDev_old {
	my (@datavalues, $n_total, $avg, $result, $i);

	$avg = $result = 0;
	@datavalues = split /\s+/, $_[0];
	$n_total = $#datavalues;
	if ($n_total > -1) {
		$avg = 0.0;
		$result = 0.0;
		
		foreach $i (@datavalues) {
			$avg += $i;
		}
		
		$avg = $avg/($n_total + 1);
		
		foreach (@datavalues) {
			$result += ($_ - $avg) **2;
		}
		$n_total = 1 if ($n_total == 0);
		$result = sqrt($result/($n_total - 1));
	}
	return ($avg, $result, ($avg * $n_total));
}

sub GetInfoInput {
	my ($input_string) = $_[0];
	$input_string =~ s/\=//g;
	my (%rec, @valid_data, $counter, $index);
	my (@data_array) = split /\s+/, $input_string;
	
	@valid_data = ("TIME", "TEMP", "PRESS");
																												   
	if ($#data_array > -1) {
		for $counter (0 ..$#data_array) {
			for $index (@valid_data) {
				if ($data_array[$counter] =~ /$index/) {
					if ($data_array[$counter + 1] && $data_array[$counter + 1] =~ /(\d+\.\d+)/) {
						$rec{$index} = $1;
					}
				}
			}
		}
	}
																												   
	return \%rec;
}
																												   

sub ParseSanderInfoFile {
	my ($infile) = $_[0];
	
	my ($rec, %OutData, $valid_info);
	open INFILE, $infile or die "Cannot open $infile: $!\n";
	while (<INFILE>) {
		if ($_ =~ /R M S  F L U C T U A T I O N S/) {
			last;
		}elsif ($_ =~ /NSTEP\s+\=\s*(\d+)\s+(.+)/) {
			$rec = GetInfoInput($2);
			if ($rec) {
				$OutData{$1} = $rec;
				$valid_info = 1;
			}
		}
	}
	close INFILE; 
	
	die "No valid data found in $infile\n"
		if (! $valid_info);

	return \%OutData;
}

sub GetMPSimEnergy {
	my ($curr_fle, $total_bases) = @_;
	my ($data_start, $in_data, $unit_nm, $counter, $eng);
	my ($curr_val, $curr_total, %Eng, @eng_header);

	if (open INPFILE, $curr_fle) {
		while (<INPFILE>) {
			chomp;
			$in_data = $_;
			if ($in_data =~ /^\s*[Atom|GROUP]/ && $#eng_header <= 0) {
				while ($in_data =~ /\s+(\w+)/g) {
					if ( ($1 !~ /GROUP/i) && ($1 !~ /TOTAL/i) ) {
						push @eng_header, $1;
					}
				}
			} elsif ($in_data =~ /^\s*(\d+)\s+/ && $#eng_header > 0) {
				$unit_nm = $1;

				if ($unit_nm == $total_bases +1) {
					last;
				} else {
					$counter = 0;
					$curr_val = $curr_total = 0;
					while ($in_data =~ /(\-?\d+\.\d+)\s*/g && $counter <= $#eng_header) {
						$eng = $eng_header[$counter];
						$curr_val = $1;
						$curr_total += $curr_val;
						$Eng{$unit_nm}->{$eng} = $curr_val; 
				$counter++;
					}								
					$Eng{$unit_nm}->{"TOTAL"} = $curr_total;
					$curr_val = $curr_total = 0;
					} 
				}
			}
		close INPFILE;
		}
		
	die "Error when reading energy file: $curr_fle: $!\n"
		if (! ($#eng_header > 0));
	
	return \%Eng;
}

sub LoadElements {
	my (%ELEMENTS, $indata, $eleNum, $myDir, $tmp, $nbnd, $nshell);
	$myDir = "$Bin/Packages";
	my ($datFile) = "$myDir/elementList.txt";
	FileTester($datFile);

	open INDATA, $datFile or die "Cannot open file $datFile: $!\n";
	while (<INDATA>) {
		chomp;
		$indata = $_;
		if ($indata =~ /^(\d+)\s+(\*?)\s*\d+/) {
			$eleNum = $1;
			if ($2) {
				$ELEMENTS{$1}{"NATURAL"} = 0;
			} else {
				$ELEMENTS{$1}{"NATURAL"} = 1;
			}
			#$indata = $'; 
			if ($indata =~ /^\d+\s+\*?\s*(\d+\.?\d*)\s+(\w+)\s+(\w+)\s+(.+)$/) {
				$ELEMENTS{$eleNum} = (
										{
											"MASS"      => $1,
											"NAME"      => $2,
											"SYMBOL"    => $3,
											"NVALENCE"  => 0,
											"NBONDS"    => 0,
											"LONEPAIRS" => 0,
										}
									);
				$tmp = $4;
				$nbnd = 0;
				while($tmp =~ /\d+(s|p|d|f)(\d+)/g) {
					if ($1 eq "s" ) {
						$nshell = 2;
					} elsif ($1 eq "p") {
						$nshell = 6;
					} elsif ($1 eq "d") {
						$nshell = 10;
					} elsif ($1 eq "f") {
						$nshell = 14;
					}
					$nbnd = $nshell - $2;
					$ELEMENTS{$eleNum}{NVALENCE} += $2;
					$ELEMENTS{$eleNum}{NBONDS} += $nbnd;
					$ELEMENTS{$eleNum}{LONEPAIRS} += int(($nshell-2*$nbnd)/2);
				}
			} else {
				delete $ELEMENTS{$eleNum};
			}
		}
	}

	close INDATA;

	die "ERORR: No valid data found in file $datFile\n"
		if (! %ELEMENTS);

	return \%ELEMENTS;

}

sub val {
	my ($inStr) = $_[0];
	my ($returnStr);

	if (IsInteger($inStr) or IsDecimal($inStr)) {
		$returnStr = $inStr;
	} elsif ($inStr =~ /^(-?\d+\.*\d*)E(\d+)/) {
		$returnStr = $1 * 10**$2;
	} elsif ($inStr =~ /^(-?\d+\.*\d*)E(\+\d+|\-\d+)/) {
		$returnStr = $1 * 10**$2;
	} else {
		$returnStr = $inStr;
	}

	return $returnStr;
}

sub CombineMols {
	my ($mol1, $mol2, $CONN, $connections, $updateRes) = @_;
	my ($atom, $tot_atoms, $tot_res, @tmp, $ATOMS, @BONDS, %CONS, $bond);

	$updateRes = 1 if (! defined($updateRes));
	@tmp = sort numerically keys %{ $mol1 };
	$tot_atoms = $tmp[$#tmp];
	$tot_res = $mol1->{$tot_atoms}{"RESNUM"};
	@tmp = sort numerically keys %{ $mol2 };

	$ATOMS = $mol1;
	%CONS = %{ $CONN };
	for $atom (@tmp) {
		$ATOMS->{($atom + $tot_atoms)} = dclone($mol2->{$atom});
		$ATOMS->{($atom + $tot_atoms)}{"RESNUM"} += $tot_res if($updateRes);
		$ATOMS->{($atom + $tot_atoms)}{"INDEX"} = $atom + $tot_atoms;
		@BONDS = ();
		@BONDS = @{ $connections->{$atom} } if (defined($connections->{$atom}));
		for $bond (@BONDS) {
			push @{ $CONS{($atom + $tot_atoms)} }, ($bond + $tot_atoms);
		}
	}

	return ($ATOMS, \%CONS);
}

sub GetBondLength {
	my ($a1, $a2, $box, $show_disp) = @_;
	my ($r, $i);

	if (! defined($box)) {
		for $i ("X", "Y", "Z") {
			$r += ($a1->{"${i}COORD"}-$a2->{"${i}COORD"})**2;
		}
		return sqrt($r);
	}

	&GetMinDist($a1,$a2,$box,\%{ $r });
	for $i ("X", "Y", "Z") {
		$r->{r2} += $r->{$i}**2;
	}

	$r->{r} = sqrt($r->{r2});

	if (! $show_disp) {
		return $r->{r};
	}else {
		return ($r->{r},$r->{disp});
	}
}

sub GetBondLengthold(@) {
	my ($a, $b, $box, $show_disp) = @_;
	my ($dx,$dy,$dz,$dx2,$dy2,$dz2,$r2);
	my ($lx,$ly,$lz,$dx2p,$dx2n);
	my ($dy2p,$dy2n,$dz2p,$dz2n);
	my ($disp, $i);

	$show_disp = 0 if (! defined($show_disp));
	$dx = $a->{XCOORD}-$b->{XCOORD};
	$dy = $a->{YCOORD}-$b->{YCOORD};
	$dz = $a->{ZCOORD}-$b->{ZCOORD};
	$dx2 = $dx*$dx;
	$dy2 = $dy*$dy;
	$dz2 = $dz*$dz;
		$disp = ();
	if(defined($box)) {
		for ("DISPX","DISPY","DISPZ") {
			$disp->{$_} = 0;
		}
		$lx = $box->{X}{len};
		$ly = $box->{Y}{len};
		$lz = $box->{Y}{len};

		$dx = fmod($dx,$lx);
		$dy = fmod($dy,$ly);
		$dz = fmod($dz,$lz);
		$dx2 = $dx*$dx;
		$dy2 = $dy*$dy;
		$dz2 = $dz*$dz;

		$dx2p = $dx2 + $lx*$lx + 2*$dx*$lx;
		$dx2n = $dx2 + $lx*$lx - 2*$dx*$lx;
		$dy2p = $dy2 + $ly*$ly + 2*$dy*$ly;
		$dy2n = $dy2 + $ly*$ly - 2*$dy*$ly;
		$dz2p = $dz2 + $lz*$lz + 2*$dz*$lz;
		$dz2n = $dz2 + $lz*$lz - 2*$dz*$lz;

		if($dx2p<$dx2) {
			$dx2 = $dx2p;
			$disp->{DISPX} = 1;
		}
		if($dx2n<$dx2) {
			$dx2 = $dx2n;
			$disp->{DISPX} = -1;
		}
		if($dy2p<$dy2) {
			$dy2 = $dy2p;
			$disp->{DISPY} = 1;
		}
		if($dy2n<$dy2) {
			$dy2 = $dy2n;
			$disp->{DISPY} = -1;
		}
		if($dz2p<$dz2) {
			$dz2 = $dz2p;
			$disp->{DISPZ} = 1;
		}
		if($dz2n<$dz2) {
			$dz2 = $dz2n;
			$disp->{DISPZ} = -1;
		}
	}
	$r2 = $dx2+$dy2+$dz2;

	if (! $show_disp) {
		return sqrt($r2);
	}else {
		return (sqrt($r2),$disp);
	}
}

sub GetAngle(@) {
	my ($bead_1, $bead_2, $bead_3, $box, $saveRad) = @_;
	my ($vec_1, $vec_2, $angle, $counter);

	if (defined($box)) {
		&GetMinDist($bead_1,$bead_2,$box,\%{ $vec_1 });
		&GetMinDist($bead_3,$bead_2,$box,\%{ $vec_2 });
	} else {
		for $counter ("X", "Y", "Z") {
			$vec_1->{$counter} = $bead_1->{"${counter}COORD"} - $bead_2->{"${counter}COORD"};
			$vec_2->{$counter} = $bead_3->{"${counter}COORD"} - $bead_2->{"${counter}COORD"};
		}
	}

# Now to calculate the dot product between the two vectors, us the dot product
# A.B = |A||B|cos(O)= sum(AiBi)

	$angle = DotProduct($vec_1, $vec_2);
	$angle = acos($angle);

	if (defined($saveRad) and $saveRad eq "1") {
		return $angle;
	} else {
		return CRadDegrees($angle, 1);
	}
}

sub GetTorsion(@) {
	my ($b1, $b2, $b3, $b4, $box, $saveRad) = @_;
	my ($counter, $v1, $v2, $v3, $v4, $v5);
	my ($normal_a, $normal_b, $answer, $sign);
	my ($pi) = atan2(1,1) *4;

	if (! defined($box)) {
		&GetMinDist($b3,$b2,$box,\%{ $v1 });
		&GetMinDist($b4,$b3,$box,\%{ $v2 });
		&GetMinDist($b2,$b3,$box,\%{ $v3 });
		&GetMinDist($b1,$b2,$box,\%{ $v4 });
		&GetMinDist($b3,$b4,$box,\%{ $v5 });
	} else {
		for $counter ("X", "Y", "Z") {
			$v1->{$counter} = $b3->{"${counter}COORD"} - $b2->{"${counter}COORD"};
			$v2->{$counter} = $b4->{"${counter}COORD"} - $b3->{"${counter}COORD"};
			$v3->{$counter} = $b2->{"${counter}COORD"} - $b3->{"${counter}COORD"};
			$v4->{$counter} = $b1->{"${counter}COORD"} - $b2->{"${counter}COORD"};
			$v5->{$counter} = $b3->{"${counter}COORD"} - $b4->{"${counter}COORD"};
		}
	}

	$normal_a = CrossProduct($v2, $v1);
	$normal_b = CrossProduct($v3, $v4);

	$answer = DotProduct($normal_a, $normal_b);
	$sign = DotProduct($normal_b, $v5);
	$answer = acos($answer);
	if ($sign < 0) {
		$answer = (2 * $pi) - $answer;
	}

	if (defined($saveRad) and $saveRad eq "1") {
		return $answer;
	} else {
		$answer = CRadDegrees($answer, 1);

		return $answer;
	}

}

sub CrossProduct(@) {
	my ($n1, $n2) = @_;

	my ($result);

	$result->{X} = ($n1->{Y} * $n2->{Z}) - ($n1->{Z} * $n2->{Y});
	$result->{Y} = ($n1->{Z} * $n2->{X}) - ($n1->{X} * $n2->{Z});
	$result->{Z} = ($n1->{X} * $n2->{Y}) - ($n1->{Y} * $n2->{X});

	return $result;

}

sub DotProduct(@) {
	my ($n1, $n2) = @_;
	my ($result, $counter, $norm_a, $norm_b);

	$result = 0;
	for $counter ("X", "Y", "Z") {
		$result += ($n1->{$counter} * $n2->{$counter});
		$norm_a += $n1->{$counter}**2;
		$norm_b += $n2->{$counter}**2;
	}

	$result = ($result/(sqrt($norm_a) * sqrt($norm_b))) if ($norm_a > 0 and $norm_b > 0);

	return $result;
}

sub VecLen {
	my ($v1) = $_[0];
	my ($len, $len_sq, $i);

	$len = $len_sq = 0.0;
	for $i ("X", "Y", "Z") {
		$len_sq += $v1->{$i}*$v1->{$i};
	}
	$len = sqrt($len_sq);
	return $len;
}

sub Normalize {
	my ($v_a) = @_;
	my ($len_a, $i);
	
	$len_a = VecLen($v_a);
	for $i ("X", "Y", "Z") {
		$v_a->{$i} /= $len_a;
	}
}

sub CRadDegrees(@) {

	my ($convertToDegrees) = $_[1];
	my ($inital_angle) = $_[0];
	my ($resulting_angle) = 0.0;

	my ($pi) = atan2(1,1) *4;

	if ($convertToDegrees) { $resulting_angle = $inital_angle * 180 / $pi; }
	else { $resulting_angle = $inital_angle * $pi/180; }

	return $resulting_angle;
}

sub acos(@) {
	my($x) = $_[0];
	if (abs($x) > 1.0) {
		return 0;
	} else {
		return atan2(sqrt(1 - $x * $x), $x);
	}
}

sub GetRotations {
	my ($axes, @vals) = @_;
	my ($rotVec);

	$rotVec = ([0,0,0]);
	$rotVec->[0] = pop @vals if ($axes =~ /x/i && $#vals > -1 && $vals[0] =~ /\-?\d+\.\d*/);
	$rotVec->[1] = pop @vals if ($axes =~ /y/i && $#vals > -1 && $vals[0] =~ /\-?\d+\.\d*/);
	$rotVec->[2] = pop @vals if ($axes =~ /z/i && $#vals > -1 && $vals[0] =~ /\-?\d+\.\d*/);
	return $rotVec;
};

sub GetRotationMatrix {
	my ($rotationVector) = @_;
	my ($count, $rotationMatrix, @cosA, @sinA);

	$count = 0;
	foreach $count (0 .. $#{ $rotationVector }) {
		$sinA[$count] = sin($rotationVector->[$count]);
		$cosA[$count] = cos($rotationVector->[$count]);
	}

   $rotationMatrix = [
						 [$cosA[1]*$cosA[2], $sinA[0]*$sinA[1]*$cosA[2]+$cosA[0]*$sinA[2], -$cosA[0]*$sinA[1]*$cosA[2]+$sinA[0]*$sinA[2]],
						 [-$cosA[1]*$sinA[2], -$sinA[0]*$sinA[1]*$sinA[2]+$cosA[0]*$cosA[2], $cosA[0]*$sinA[1]*$sinA[2]+$sinA[0]*$cosA[2]],
						 [$sinA[1], -$sinA[0]*$cosA[1], $cosA[0]*$cosA[1]],
					 ];
	return $rotationMatrix;
};


sub Rotate {
   my ($atoms, $rot_array, $coord) = @_;
   my (@rotation_matrix, $counter, $index, $matrix_keys, %F_Data, $i, $offset);
   my (@sinA, @cosA);

   for $counter (keys %{ $atoms }) {
	   for $index (keys %{ $atoms->{$counter} }) {
		   $F_Data{$counter}{$index} = $atoms->{$counter}{$index};
	   }
   }

   for $counter (0 .. $#{ $rot_array }) {
		$sinA[$counter] = sin($rot_array->[$counter]);
		$cosA[$counter] = cos($rot_array->[$counter]);
   }

   @rotation_matrix = (
					   {
						   "XCOORD" => [1, 0, 0],
						   "YCOORD" => [0, $cosA[0], $sinA[0]],
						   "ZCOORD" => [0, -$sinA[0], $cosA[0]],
						},

					   {
						   "XCOORD" => [$cosA[1], 0, -$sinA[1]],
						   "YCOORD" => [0, 1, 0],
						   "ZCOORD" => [$sinA[1], 0, $cosA[1]],
					   },
					   {
						   "XCOORD" => [$cosA[2], $sinA[2], 0],
						   "YCOORD" => [-$sinA[2], $cosA[2], 0],
						   "ZCOORD" => [0, 0, 1],
					   },
					   {
						   "XCOORD" => [$cosA[1]*$cosA[2], $sinA[0]*$sinA[1]*$cosA[2]+$cosA[0]*$sinA[2], -$cosA[0]*$sinA[1]*$cosA[2]+$sinA[0]*$sinA[2]],
						   "YCOORD" => [-$cosA[1]*$sinA[2], -$sinA[0]*$sinA[1]*$sinA[2]+$cosA[0]*$cosA[2], $cosA[0]*$sinA[1]*$sinA[2]+$sinA[0]*$cosA[2]],
						   "ZCOORD" => [$sinA[1], -$sinA[0]*$cosA[1], $cosA[0]*$cosA[1]],
					   },
					   );

   for $matrix_keys (keys %{ $atoms }) {
	   $counter = \%{ $rotation_matrix[$coord] };
	   for $index ("XCOORD", "YCOORD", "ZCOORD") {
		   $atoms->{$matrix_keys}{$index} = ($F_Data{$matrix_keys}{"XCOORD"} * $counter->{$index}[0] +
											$F_Data{$matrix_keys}{"YCOORD"} * $counter->{$index}[1] +
											$F_Data{$matrix_keys}{"ZCOORD"} * $counter->{$index}[2]);
	   }
   }
   @rotation_matrix = ();

}

sub CenterMol {
   my ($atoms) = @_;
   my ($CENTER, $atomC, $dim);

	$CENTER = CoP($atoms);
	for $atomC (keys %{ $atoms }) {
		for $dim ("XCOORD", "YCOORD", "ZCOORD") {
			$atoms->{$atomC}{$dim} -= $CENTER->{$dim};
		}
   }

   return $CENTER;
}

sub CoP {
	my ($atoms) = @_;
	my (%CENTER, $atomC, $tot, $dim);
	
   $tot = 0;
   for $atomC (keys %{ $atoms }) {
		for $dim ("XCOORD", "YCOORD", "ZCOORD") {
			$CENTER{$dim} += $atoms->{$atomC}{$dim};
		}
		$tot++;
   }

   for $dim ("XCOORD", "YCOORD", "ZCOORD") {
		$CENTER{$dim} /= $tot;
   }

	return \%CENTER;
}

sub PrintProgress {
	my ($currPos, $total, $start, $pStr) = @_;
	my ($progress, $str, $end);
	
	$end = time();

	$progress = $currPos/$total;
	
	$str = sprintf("%.2f%% complete %s\r", 
				   100*$progress, getEta(($end - $start), $progress));
	
	print "${pStr}${str}" if (defined($pStr));
	return length($str);
}

sub GetTime {
	my ($timeLeft) = $_[0];
	my ($returnStr);
 
	if ($timeLeft > 60) {
		if ($timeLeft > 3600) {
			if ($timeLeft > 86400) {
				$returnStr = int($timeLeft/86400) . "d ";
				$timeLeft = $timeLeft % 86400;
			}
			$returnStr .= int($timeLeft/3600) . "h ";
			$timeLeft = $timeLeft % 3600;
		}
		$returnStr .= int($timeLeft/60) . "m ";
		$timeLeft = $timeLeft % 60;
	}
	$returnStr .= sprintf("%.0f", $timeLeft);
	return $returnStr;
}

sub getEta {
	my ($elapsed, $percentage) = @_;
	my ($totalTime) = $elapsed/$percentage;
	my ($timeLeft) = $totalTime - $elapsed;
	my ($returnStr) = "(";
	
	$returnStr .= GetTime($timeLeft) . "s remaining)		  ";
	
	return $returnStr;
}

sub GetSelections {
	my ($opts, $isPrint) = @_;
	my ($index, %SELECT, $saveKey, $counter, $exclude, %EXCLUDED, $interval, @tmp, $i);

	$isPrint = 1 if (! defined($isPrint));

	for $index (@{ $opts }) {
		$index = Trim($index);
		$exclude = 0;
		if ($index =~ /^\^/) {
			$exclude = 1;
			$index =~ s/^\^//;
		}
		if ($index eq "*" or lc($index) eq "all") {
			$SELECT{INDEX}{"*"} = 1;
		} elsif ($index =~ /^([A-Z])([A-Z])(\S+)\:?/i) {
			$saveKey = getOptName($1, $2);
			if (defined($saveKey)) {
				print "SELECTED $saveKey $3\n" if ($isPrint);
				$counter = $3;
				$counter = "!" . $counter if ($exclude);
				$SELECT{$saveKey}{$counter} = 1;
			}
		} elsif ($index =~ /^\:([A-Z])([A-Z])(\d+)\-(\d+)\:?(\d*)(r)?/i) {
			$saveKey = getOptName($1, $2);
			$interval = 1;
			$interval = $5 if ($5 ne "");
			if (defined($saveKey)) {
				if (! $6) {
					print "SELECTED $saveKey $3 - $4 every $interval\n" if ($isPrint);
					$counter = $3;
					while ($counter <= $4) {
						$SELECT{$saveKey}{$counter} = 1;
						$counter += $interval;
					}
				} else {
					print "Will select $interval random $saveKey between $3 - $4\n" if ($isPrint);
					for ($3 .. $4) {
						push @tmp, $_;
					}
					$counter = 1;
					while ($counter <= $interval or $#tmp == -1) {
						$i = int(rand($#tmp));
						$SELECT{$saveKey}{$tmp[$i]} = 1;
						splice @tmp, $i, 1;
						$counter++;
					}
				}
			}
		}
	}

	die "ERROR: No valid options found!\n"
		if (! %SELECT);

	return \%SELECT;
}

sub getOptName {
	my ($optField, $optType) = @_;
	my ($saveKey);

	($optField, $optType) = (uc $optField, uc $optType);
	
	if ($optField eq "C" and $optType =~ /(X|Y|Z)/i) {
		$saveKey = "${optType}COORD";
		$saveKey = uc($saveKey);
	}elsif (uc($optField) eq "N") {
		if (uc($optType) eq "A") {
			$saveKey = "ATMNAME";
		} elsif (uc($optType) eq "R") {
			$saveKey = "RESNAME";
		} elsif (uc($optType) eq "M") {
			$saveKey = "MOLECULEID";
		}
	}elsif ($optField eq "T" and $optType eq "A") {
		$saveKey = "FFTYPE";
	}elsif ($optField eq "T" and $optType eq "C") {
		$saveKey = "CHAIN";
	}elsif ($optField eq "S" and $optType eq "M") {
		$saveKey = "MOLSIZE";
	} else {
		if ($optType eq "A") {
			$saveKey = "INDEX";
		} elsif ($optType eq "R") {
			$saveKey = "RESNUM";
		} elsif ($optType eq "M") {
			$saveKey = "MOLECULEID"
		} else {
			$saveKey = "NUMBONDS";
		}
	}
	
	return $saveKey;
}

sub printAtomSelectOptions {
	print <<"END";
		[^][:][I|T|N|S|C][a|r] [>|<]
			a   - atom
				r   - residue
				IaX - atom number X
				IrX - residue index X
				TaX - atom type X
				NaX - atom name X
				NrX - residue name X
				SmX - molecule size
				CXYZ -  atom X|Y|Z coordinate
				Use ":" to specify a range, eg. :Tr1-8 :Ia3-66
				Use "^" to exclude a selection. You can use multiple combinations
				range and exclusion enclosed in quotes, eg, "^:TrIP-IM ^:Ia23-45"
				to exclude residues of type IM and IP and atoms 23-45
				Use ">" or "<" for less than and greater than respectively
END
}
	
sub CoM {
	my ($atoms) = $_[0];
	my ($i, %COM, $dim, $totalMass, $mass, @tmp, $ori, $F);
	
	$totalMass = 0;
	@tmp = keys %{ $atoms };
	return $atoms->{$tmp[0]} if ($#tmp == 0);
	for $i (@tmp) {
		$mass = 1;
		$mass = $atoms->{$i}{"MASS"} if (exists($atoms->{$i}{MASS}));
		for $dim ("XCOORD", "YCOORD", "ZCOORD") {
			$COM{$dim} += $atoms->{$i}{$dim} * $mass;
		}
		$totalMass += $mass;
	}
	
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
		$COM{$dim} /= $totalMass;
	}
   
	return \%COM;
}

sub ParseSelect {
	my ($selection) = $_[0];
	my ($i, %SELECT);

	print "SNAPSHOT SELECTION: ";
	if (defined($selection) && Trim($selection) =~ /^(\d+)\s+(\d+)\s+(\d+)/) {
		print "start: $1 end: $2 every: $3...";
		$i = $1;
		while ($i <= $2) {
			$SELECT{$i} = 1;
			$i += $3;
		}
	} else {
		%SELECT = ();
		print "using all";
	}
	print "\n";

	return \%SELECT;
}

sub execCmd {
	my ($cmd, $fileToDelete) = @_;

	die "Error while executing $cmd\n" if (system($cmd));
	system "rm -f $fileToDelete" if (defined($fileToDelete));
}

sub TrjSelections {
	my ($select) = $_[0];
	my (@tmp, $i, %SELECTION, $j, $counter, @tmp2);

	@tmp = split /\s+/, $select;

	$i = 0;
	$counter = 0;
	while ($i <= $#tmp) {
		if ($tmp[$i] =~ /^\:It(\d+)\-(\d+)\:(\d+)$/) {
			$j = $1;
			while ($j <= $2) {
				$SELECTION{$j} = 1;
				$counter++;
					$j += $3;
			}
		} elsif ($tmp[$i] =~ /^(\d+)$/) {
			$SELECTION{$1} = 1;
			$counter++;
		} elsif (-e $tmp[$i] && -r $tmp[$i] && -T $tmp[$i]) { # if the input is a file
			if (open INDATA, $tmp[$i]) {
				while (<INDATA>) {
					if ($_ =~ /^\"(.+)\"/) {
						@tmp2 = split /\s+/, $1;
						for $j (@tmp2) {
							push @tmp, $j;
						}
					}
				}
				close INDATA;
			}
		} elsif ($tmp[$i] =~ /^last$/i) {
			$SELECTION{"-1"} = 1;
			$counter++;
		}
		$i++;
	}

	die "ERROR: No valid selection found. See help file. Got $select\n" if ($select ne "*" && (! %SELECTION || ! keys %SELECTION));

	print "trajectory: ";
	if ($select eq "*") {
		print "selected all frames...";
	} else {
		print "selected $counter frames...";
	}
	return \%SELECTION;
}

sub CenterOnMol {
	my ($atoms, $CENTER) = @_;
	my ($i, @tmp, $j);

	@tmp = ("XCOORD", "YCOORD", "ZCOORD");
	for $i (keys %{ $atoms }) {
		for $j (@tmp) {
			$atoms->{$i}{$j} -= $CENTER->{$j};
		}
	}
}

sub GetSoluteAtoms {
	my ($atoms, $moleculeList) = @_;
	my (%SOL, $i, $atm, @tmp);

	for $i (keys %{ $moleculeList }) {
		@tmp = keys %{ $moleculeList->{$i}{MEMBERS} };
		$atm = pop @tmp;
		next if ($atoms->{$atm}{RESNAME} =~ /WAT|Na|Mg|Cl/i || $atoms->{$atm}{FFTYPE} =~ /OF3C|Cl|Na|Mg|HF3C|OW|HW/i);
		$SOL{$atm} = 1;
		for $atm (@tmp) {
			$SOL{$atm} = 1;
		}
	}

	return \%SOL;
}

sub Permutate {
	my @head = @{ $_[0] };
	my @tail = @{ $_[1] };
	my (@RET);
	
	@RET = ();
	unless (@head) {
		# stop recursing when there are no elements in the head
		return \@tail;
	} else {
		# for all elements in @head, move one from @head to @tail
		# and call permut() on the new @head and @tail
		my(@newhead,@newtail,$i);
		foreach $i (0 .. $#head) {
			@newhead = @head;
			@newtail = @tail;
			unshift(@newtail, splice(@newhead, $i, 1));
			push @RET, Permutate([@newhead], [@newtail]);
		}
	}

	return @RET;
}

sub ShowSelectionInfo {
	my ($usage) = <<DATA;
Selection options:
		[^][:][I|T|N][a|r|b|m]
		a   - atom
		r   - residue
		m   - # atoms in molecule (only valid with index)
		b   - total number of bonds of atom (only valid with index)
		IaX - atom number X
		IrX - residue index X
		TaX - atom type X
		NaX - atom name X
		NrX - residue name X
		Use ":" to specify a range, eg. :Tr1-8 :Ia3-66
		Use "^" to exclude a selection. You can use multiple combinations
		range and exclusion enclosed in quotes, eg, "^:TrIP-IM ^:Ia23-45"
		to exclude residues of type IM and IP and atoms 23-45
DATA

	return $usage;
}

sub ReadParmFile {
	my ($parmFile, $PARMS) = @_;
	my ($valid, $i);

	$valid = 0;
	open PARMFILE, $parmFile || die "ERROR: Cannot access param file $parmFile: $!\n";
	while (<PARMFILE>) {
		chomp;
		if ($_ =~ /^PARM_(\w+)\s+(.+)$/i) {
			$PARMS->{uc $1}{VAL} = $2;
			$valid = 1;
		}
	}
	close PARMFILE;

	die "ERROR: No valid data found while reading parm file $parmFile!\n" if (! $valid);

	for $i (keys %{ $PARMS }) {
		die "ERROR: $i is required but not found in parm file!\n" 
			if (exists($PARMS->{$i}{REQUIRED}) and ! exists($PARMS->{$i}{VAL}));
	}
}

sub ShuffleArray {
# fisher_yates_shuffle
	my $array = $_[0];
	my $i;
	for ($i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}
}

sub GetMinDist {
	my ($a1, $a2, $box, $pv) = @_;
	my ($i, @r, @ori, @cb,@n, $disp);

	$i = 0;
	$disp = ();
	for ("X","Y","Z") {
		$r[$i] = $a1->{"${_}COORD"}-$a2->{"${_}COORD"};
		$cb[$i] = $box->{$_}{lo} + $box->{$_}{len}/2;
		$i++;
	}
	$ori[0]=$r[0]+$cb[0]; #move to the center of box
	$ori[1]=$r[1]+$cb[1]; 
	$ori[2]=$r[2]+$cb[2];
	$n[0]=$box->{Hinv}[0][0]*$ori[0]+$box->{Hinv}[0][1]*$ori[1]+$box->{Hinv}[0][2]*$ori[2];
	$n[1]=$box->{Hinv}[1][0]*$ori[0]+$box->{Hinv}[1][1]*$ori[1]+$box->{Hinv}[1][2]*$ori[2];
	$n[2]=$box->{Hinv}[2][0]*$ori[0]+$box->{Hinv}[2][1]*$ori[1]+$box->{Hinv}[2][2]*$ori[2];

	while($n[0]>1) { $n[0] -= 1; $disp->{DISPX} += 1; }
	while($n[0]<0) { $n[0] += 1; $disp->{DISPX} -= 1; }
	while($n[1]>1) { $n[1] -= 1; $disp->{DISPY} += 1; }
	while($n[1]<0) { $n[1] += 1; $disp->{DISPY} -= 1; }
	while($n[2]>1) { $n[2] -= 1; $disp->{DISPZ} += 1; }
	while($n[2]<0) { $n[2] += 1; $disp->{DISPZ} -= 1; }

	$pv->{X}=$box->{H}[0][0]*$n[0]+$box->{H}[0][1]*$n[1]+$box->{H}[0][2]*$n[2]-$cb[0];
	$pv->{Y}=$box->{H}[1][0]*$n[0]+$box->{H}[1][1]*$n[1]+$box->{H}[1][2]*$n[2]-$cb[1];
	$pv->{Z}=$box->{H}[2][0]*$n[0]+$box->{H}[2][1]*$n[1]+$box->{H}[2][2]*$n[2]-$cb[2];
	$pv->{disp} = $disp;
}

sub dPcmplx {
	my ($v1, $v2, $ccfactor) = @_;
	my ($i, $val);

	$val = 0;
	for $i (keys %{ $v1 }) {
		$val += $v1->{$i}{re}*$v2->{$i}{re} + $ccfactor*$v1->{$i}{im}*$v2->{$i}{im}; 
		#note that dotP(v1,v2) = sum(v1i*cc(v2i))
		#where cc is complex conjugate
		#so there are 2 sign flips that cancel
		#i^2 and -1*im(b)
	}
	return $val;
}

sub HasCell {
	my ($headers) = $_[0];
	my ($hasCell, $i);

	$hasCell = 0;

	for $i (@{ $headers }) {
		if ($i =~ /^CRYSTX\s+(.+)/) {
			$hasCell = 1;
			last;
		}
	}

	return $hasCell;
}

sub AddElementField {
	my ($atoms, $elements) = @_;
	my ($i, $sfield, $elementNum);

	for $i (keys %{ $atoms }) {
		$sfield = $atoms->{$i}{FFTYPE};
		$sfield =~ s/_.*//;
		$sfield =~ s/\d+.*//;
		$elementNum = FindElement($elements, $sfield);
		if (defined($elementNum)) {
			$atoms->{$i}{ELEMENT} = \%{ $elements->{$elementNum} };
		} else {
			$sfield = $atoms->{$i}{ATMNAME};
			$sfield =~ s/_.*//;
			$sfield =~ s/\d+.*//;
			$elementNum = FindElement($elements, $sfield);
			if (! defined($elementNum)) {
				die "ERROR: Cannot find element for atom # $i: $sfield\n";
			}
			$atoms->{$i}{ELEMENT} = \%{ $elements->{$elementNum} };
		}
	}
}

sub FindElement {
	my ($elements, $sfield) = @_;
	my ($i, $eleNum,$sym);

	for $i (keys %{ $elements }) {
		$sym = $elements->{$i}{SYMBOL};
		if ($sfield =~ /^${sym}/i) {
			$eleNum = $i;
			last;
		}
	}
	die "ERROR: Cannot find element type of fftype $sfield\n"
		if (!defined($eleNum));
	return $eleNum;
}

sub AddRingField {
	my ($atoms, $bonds, $select) = @_;
	my ($i, $ring, $cycle, $path);

	for $i (keys %{ $select }) {
		$cycle = $path = ();
		&getMinCycle($bonds, \@{ $cycle }, \@{ $path }, $i, $i);
		$atoms->{$i}{RING} = 0 if ($#{ $cycle } == -1);
		$atoms->{$i}{RING} = $#{ $cycle } + 1 if ($#{ $cycle } > -1);
		@{ $atoms->{$i}{CYCLE} } = @{ $cycle };
	}
}

sub getMinCycle {
	my ($bonds, $minCycle, $currPath, $i, $j) = @_;
	my ($k, $l, $parent);

	$parent = $currPath->[ $#{ $currPath }]; #last element of path array is parent
	push @{ $currPath }, $i; 
	for $k (@{ $bonds->{$i} }) {
		next if ($k == $parent); #prevent trivial backtracking
		if($k == $j) { #cycle
			if ($#{ $minCycle } == -1 or $#{ $currPath } < $#{ $minCycle }) { #min cycle check
				@{ $minCycle } = @{ $currPath }; #found smaller cycle, so assign and return to previous
				pop @{ $currPath };
				return;
			}
		} elsif (alreadyVisited($currPath,$k, $j)) { #not a cycle if already visited, try next path(?)
			next;
		} 
		#else continue along path recursively...
		&getMinCycle($bonds, $minCycle, $currPath, $k, $j);
	}
	pop @{ $currPath }; #not a cycle, so backtrack to parent and try another path(?)
}

sub alreadyVisited {
	my ($nlist, $cnode, $snode) = @_;
	my ($visited, $i, $tmp);

	%{ $tmp } = map { $_ => 1 } @{ $nlist };
	return 1 if(exists($tmp->{$cnode}));
	return 0;
}

1;
