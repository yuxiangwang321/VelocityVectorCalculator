proc vvc {file_1 file_2 outputName baseFrame step type xl xh yl yh zl zh binSize scale atomSelection} {

# output velocity according displacement between two frames, step cannot be 0
# type can be: x y z xy xz yz xyz
# scale is the time between two continous frames in the original trajectory file.
# velocity unit is: nm/ns
# angle unit in output file is degree.
# example command: source VelocityVectorCalculator.tcl; vvc ce0-wt.psf main.dcd water_flow 0 10 xy 0 120 0 110 0 120 3.0 0.005 {name OT}

# shared code------------------------------------------------------------
	set targetFrame [expr {$baseFrame+$step}]
	if {[string equal $file_1 $file_2]} {
		mol new $file_1 first $baseFrame last $targetFrame step $step waitfor all
	} else {
		mol new $file_1
		mol addfile $file_2 first $baseFrame last $targetFrame step $step waitfor all
	}

	set tcl_precision 8
	set pi 3.14159
	# pbc unwrap -all

	# slice the cell into bins-------------------------------------------
	set binNumReal [expr double($xh-$xl)/$binSize]
	set binNumInteger [expr int(($xh-$xl)/$binSize)]
	if { $binNumReal > $binNumInteger } {
		set binNumX [expr {$binNumInteger+1}]
	} else {
		set binNumX $binNumInteger
	}

	set binNumReal [expr double($yh-$yl)/$binSize]
	set binNumInteger [expr int(($yh-$yl)/$binSize)]
	if { $binNumReal > $binNumInteger } {
		set binNumY [expr {$binNumInteger+1}]
	} else {
		set binNumY $binNumInteger
	}

	set binNumReal [expr double($zh-$zl)/$binSize]
	set binNumInteger [expr int(($zh-$zl)/$binSize)]
	if { $binNumReal > $binNumInteger } {
		set binNumZ [expr {$binNumInteger+1}]
	} else {
		set binNumZ $binNumInteger
	}
	# slice the cell into bins --- END
# shared code --- END

# 1D velocity distribution in Y axis----------------------------------------------
	if {[string equal -nocase $type y]} {
		set filename $outputName
		set writefile [open $filename.dat w+]
		puts $writefile "#------Velocity distribution along $type axis ------"
		puts $writefile "#Command: v123 $file_1 $file_2 $baseFrame $step $type $xl $xh $yl $yh $zl $zh $binSize $scale {$atomSelection}"
		puts $writefile "#	x_init y_init x_end y_end vx vy v_mag"
		puts "------Velocity distribution along $type axis ------"
		puts "#	x_init y_init x_end y_end vx vy v_mag"
		
		for {set j 0} {$j<$binNumY} {incr j 1} {
			set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
					$xl $xh [expr {$yl+$binSize*$j}] [expr {$yl+$binSize*($j+1)}]  $zl $zh]			
			set binAtomInit [atomselect top "$bin and $atomSelection" frame 0]
			set binAtomIndex [$binAtomInit get index]
			set binLength [llength $binAtomIndex]
			
			if {$binLength > 0} {
				set binPosInit [measure center $binAtomInit]
				set binAtomEnd [atomselect top "index $binAtomIndex" frame 1]
				set binPosEnd [measure center $binAtomEnd]
				
				set vx [expr ([lindex $binPosEnd 0]-[lindex $binPosInit 0])/$scale/$step]
				set vy [expr ([lindex $binPosEnd 1]-[lindex $binPosInit 1])/$scale/$step]
				set vz [expr ([lindex $binPosEnd 2]-[lindex $binPosInit 2])/$scale/$step]			
			} else {
				set vx 0.0
				set vy 0.0
				set vz 0.0
			}
			
			set yPos [expr {$yl+$binSize*($j+0.5)}]
			puts $writefile [format "%.2f %.4f %.4f %.4f" $yPos $vx $vy $vz]
			puts [format "%.2f %.4f %.4f %.4f" $yPos $vx $vy $vz]
		}
		close $writefile
# 1D velocity distribution in Y axis --- END

# 2D velocity distribution on XY plane----------------------------------------------
	}	elseif {[string equal -nocase $type xy]} {
		set filename $outputName
		set writefile [open $filename.dat w+]
		puts $writefile "#------Velocity distribution on $type plan ------"
		puts $writefile "#Command: v123 $file_1 $file_2 $baseFrame $step $type $xl $xh $yl $yh $zl $zh $binSize $scale {$atomSelection}"
		puts $writefile "# x_init   y_init       dx       dy    d_mag         vx         vy      v_mag"
		puts "------Velocity distribution on $type plane ------"
		puts "# x_init   y_init       dx       dy    d_mag         vx         vy      v_mag"
		
		for {set i 0} {$i<$binNumX} {incr i 1} {
			for {set j 0} {$j<$binNumY} {incr j 1} {
				set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$binSize*$i}] [expr {$xl+$binSize*($i+1)}] \
						[expr {$yl+$binSize*$j}] [expr {$yl+$binSize*($j+1)}] $zl $zh]			
				set binAtomInit [atomselect top "$bin and $atomSelection" frame 0]
				set binAtomIndex [$binAtomInit get index]
				set binLength [llength $binAtomIndex]
				
				if {$binLength > 0} {
					set binPosInit [measure center $binAtomInit]
					set binAtomEnd [atomselect top "index $binAtomIndex" frame 1]
					set binPosEnd [measure center $binAtomEnd]
					
					set x_init [lindex $binPosInit 0]
					set y_init [lindex $binPosInit 1]
					set x_end  [lindex $binPosEnd 0]
					set y_end  [lindex $binPosEnd 1]
					set dx [expr {$x_end-$x_init}]
					set dy [expr {$y_end-$y_init}]
					set d_mag [expr {sqrt($dx*$dx + $dy*$dy)}]
					
					set vx [expr ([lindex $binPosEnd 0]-[lindex $binPosInit 0])/10/$scale/$step]
					set vy [expr ([lindex $binPosEnd 1]-[lindex $binPosInit 1])/10/$scale/$step]	
					set v_mag  [expr {sqrt($vx*$vx + $vy*$vy)}]
					# set v_angle [expr {180*atan($vy/$vx)/$pi}]		

					puts $writefile [format "%8.2f %8.2f %8.2f %8.2f %8.2f %10.4f %10.4f %10.4f " $x_init $y_init $dx $dy $d_mag $vx $vy $v_mag]
					puts [format "%8.2f %8.2f %8.2f %8.2f %8.2f %10.4f %10.4f %10.4f " $x_init $y_init $dx $dy $d_mag $vx $vy $v_mag]
					
				} else {
					# do nothing
				}				
			}
		}
		close $writefile
# 2D velocity distribution on XY plane --- END

# 2D velocity distribution on YZ plane----------------------------------------------
	}	elseif {[string equal -nocase $type yz]} {
		set filename $outputName
		set writefile [open $filename.dat w+]
		puts $writefile "#------Velocity distribution on $type plan ------"
		puts $writefile "#Command: v123 $file_1 $file_2 $baseFrame $step $type $xl $xh $yl $yh $zl $zh $binSize $scale {$atomSelection}"
		puts $writefile "#	xPos	yPos	vx	vy	v_mag v_angle"
		puts "------Velocity distribution on $type plane ------"
		puts "	yPos	zPos	vx	vy	v_mag v_angle"
		
		for {set i 0} {$i<$binNumX} {incr i 1} {
			for {set j 0} {$j<$binNumY} {incr j 1} {
				set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$binSize*$i}] [expr {$xl+$binSize*($i+1)}] \
						[expr {$yl+$binSize*$j}] [expr {$yl+$binSize*($j+1)}] $zl $zh]			
				set binAtomInit [atomselect top "$bin and $atomSelection" frame 0]
				set binAtomIndex [$binAtomInit get index]
				set binLength [llength $binAtomIndex]
				
				if {$binLength > 0} {
					set binPosInit [measure center $binAtomInit]
					set binAtomEnd [atomselect top "index $binAtomIndex" frame 1]
					set binPosEnd [measure center $binAtomEnd]
					
					set x_init [lindex $binPosInit 0]
					set y_init [lindex $binPosInit 1]
					set x_end  [lindex $binPosEnd 0]
					set y_end  [lindex $binPosEnd 1]
					set dx [expr {$x_end-$x_init}]
					set dy [expr {$y_end-$y_init}]
					set d_mag [expr {sqrt($dx*$dx + $dy*$dy)}]
					
					set vx [expr ([lindex $binPosEnd 0]-[lindex $binPosInit 0])/$scale/$step]
					set vy [expr ([lindex $binPosEnd 1]-[lindex $binPosInit 1])/$scale/$step]	
					set v_mag  [expr {sqrt($vx*$vx + $vy*$vy)}]
					# set v_angle [expr {180*atan($vy/$vx)/$pi}]		

					puts $writefile [format "%8.2f %8.2f %8.2f %8.2f %8.2f %10.4f %10.4f %10.4f " $x_init $y_init $dx $dy $d_mag $vx $vy $v_mag]
					puts [format "%8.2f %8.2f %8.2f %8.2f %8.2f %10.4f %10.4f %10.4f " $x_init $y_init $dx $dy $d_mag $vx $vy $v_mag]
					
				} else {
					# do nothing
				}				
			}
		}
		close $writefile
# 2D velocity distribution on YZ plane --- END

# 3D velocity distribution on XYZ cube----------------------------------------------
	}	elseif {[string equal -nocase $type xyz]} {
		set filename $outputName
		set writefile [open $filename.dat w+]
		puts $writefile "#------Velocity distribution in $type cube ------"
		puts $writefile "#Command: v123 $file_1 $file_2 $baseFrame $step $type $xl $xh $yl $yh $zl $zh $binSize $scale {$atomSelection}"
		puts $writefile "#	yPos	xPos	vx	vy	vz	vxyz"
		puts "------Velocity distribution in $type cube ------"
		puts "yPos	xPos	vx	vy	vz	vxyz"
		
		for {set i 0} {$i<$binNumX} {incr i 1} {
			for {set j 0} {$j<$binNumY} {incr j 1} {
				for {set k 0} {$k<$binNumZ} {incr k 1} {
					set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
							[expr {$xl+$binSize*$i}] [expr {$xl+$binSize*($i+1)}] \
							[expr {$yl+$binSize*$j}] [expr {$yl+$binSize*($j+1)}] \
							[expr {$zl+$binSize*$k}] [expr {$zl+$binSize*($k+1)}]]			
					set binAtomInit [atomselect top "$bin and $atomSelection" frame 0]
					set binAtomIndex [$binAtomInit get index]
					set binLength [llength $binAtomIndex]
					
					if {$binLength > 0} {
						set binPosInit [measure center $binAtomInit]
						set binAtomEnd [atomselect top "index $binAtomIndex" frame 1]
						set binPosEnd [measure center $binAtomEnd]
						
						set vx [expr ([lindex $binPosEnd 0]-[lindex $binPosInit 0])/$scale/$step]
						set vy [expr ([lindex $binPosEnd 1]-[lindex $binPosInit 1])/$scale/$step]
						set vz [expr ([lindex $binPosEnd 2]-[lindex $binPosInit 2])/$scale/$step]			
					} else {
						set vx 0.0
						set vy 0.0
						set vz 0.0
					}
					
					set vxyz [expr sqrt($vx*$vx + $vy*$vy + $vz*$vz)]
					set xPos [expr {$xl+$binSize*($i+0.5)}]
					set yPos [expr {$yl+$binSize*($j+0.5)}]
					set zPos [expr {$zl+$binSize*($k+0.5)}]
					puts $writefile [format "%.2f %.2f %.4f %.4f %.4f %.4f" $xPos $yPos $vx $vy $vz $vxyz]
					puts [format "%.2f %.2f %.4f %.4f %.4f %.4f" $xPos $yPos $vx $vy $vz $vxyz]
				}
			}
		}
		close $writefile
	}
# 3D velocity distribution on XYZ cube --- END

} 
# end of code

