<?php

//PEAKFINDER Noise tolerant fast peak finding algorithm
//   INPUTS:
//       x0 - A real vector from the maxima will be found (required)
//       sel - The amount above surrounding data for a peak to be
//           identified (default = (max(x0)-min(x0))/4). Larger values mean
//           the algorithm is more selective in finding peaks.
//       thresh - A threshold value which peaks must be larger than to be
//           maxima or smaller than to be minima.
//       extrema - 1 if maxima are desired, -1 if minima are desired
//           (default = maxima, 1)
//   OUTPUTS:
//       peakLoc - The indicies of the identified peaks in x0
//       peakMag - The magnitude of the identified peaks
//
//   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
//       are at least 1/4 the range of the data above surrounding data.
//
//   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
//       that are at least sel above surrounding data.
//
//   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local 
//       maxima that are at least sel above surrounding data and larger
//       (smaller) than thresh if you are finding maxima (minima).
//
//   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
//       data if extrema > 0 and the minima of the data if extrema < 0
//
//   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
//       local maxima as well as the magnitudes of those maxima
//
//   If called with no output the identified maxima will be plotted along
//       with the input data.
//
//   Note: If repeated values are found the first is identified as the peak
//
// Ex:
// t = 0:.0001:10;
// x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
// x(1250:1255) = max(x);
// peakfinder(x)
//
// Copyright Nathanael C. Yoder 2011 (nyoder@gmail.com)

// Perform error checking and set defaults if not passed in
// Rewrite by Davide Gazzè 12/06/2012

function diff($x0, $min=null, $max=null) {
    $len = count($x0);
    if($min==null) $min = 1;
    if($max==null) $max = $len;
    $len = count($x0);
    $diff_array = array();
    for($i=$min; $i<$max; $i++) {
        $diff_array[$i-1] = $x0[$i] - $x0[$i-1];
    }
    return $diff_array;
}

function sign($number) {
    // Return
    //  1 : positive
    // -1 : negative
    //  0 : 0
    return ($number>0) ? 1 : (($number<0) ? -1 : 0); 
} 

function signArray($x0) {
    $len = count($x0);
    $result = array();
    for($i=0; $i<$len; $i++) {
        $result[$i] = sign($x0[$i]);
    }
    return $result;
}

function maxIndex($x) {
    $max = $x[0];
    $index = 0;
    for($i=1; $i<count($x) ; $i++) {
        if($x[$i] > $max) {
            $max = $x[$i];
            $index = $i;
        }
    }
    $response = array();
    $response['value'] = $max;
    $response['index'] = $index;
    return $response;
}

function findChangeSign(&$x0) {
    //find(dx0(1:end-1).*dx0(2:end) < 0)+1
    $lenX0 = count($x0)-1;
    $result = array();
    for($i=0; $i<$lenX0; $i++) {
        if(($x0[$i]*$x0[$i+1])<=0) {
            $result[] = $i+1;
        }
    }
    return $result;
}

function zeros($num) {
    $result = array();
    for($i=0; $i<$num; $i++) {
        $result[$i] = 0;
    }
    return $result;
}

function mulArray($x0, $mul) {
    for($i=0; $i<count($x0); $i++) {
        $x0[$i] = $x0[$i]*$mul;
    }
    return $x0;
}

function peakfindex($x0, $sel = -1, $thresh = null, $extrema = 1) {
    if($sel==-1) $sel = (max($x0)-min($x0))/4;
    // Make it so we are finding maxima regardless
    $x0 = mulArray($x0, $extrema);
    // Adjust threshold according to extrema
    $thresh = $thresh*$extrema;
    // Find derivative
    $dx0 = diff($x0); 
    // Find where the derivative changes sign
    $index_raw = findChangeSign($dx0);
    // Include endpoints in potential peaks and valleys
    // x = [x0(1);x0(ind);x0(end)];
    // ind = [1;ind;len0];
    $lenX0 = count($x0);
    $x = array();
    $index = array();
    $j=0;
    $index[$j] = 0;
    $x[$j] = $x0[0];
    $j++;
    for($i=0; $i<count($index_raw); $i++) {
        $x[$j] = $x0[$index_raw[$i]];
        $index[$j] = $index_raw[$i];
        $j++;
    }
    $x[$j] = $x0[$lenX0-1];
    $index[$j] = $lenX0-1;
    // x only has the peaks, valleys, and endpoints
    $lenX = count($x);
    $minMag = min($x);
    // Function with peaks and valleys
    $ii = 0;
    if($lenX>2) {
        // Set initial parameters for loop
        $tempMag = $minMag;
        $foundPeak = false;
        $leftMin = $minMag;
        //Deal with first point a little differently since tacked it on
        //Calculate the sign of the derivative since we taked the first point
        //on it does not neccessarily alternate like the rest.
        $signDx = signArray(diff($x, 1, 3));
        // The first point is larger or equal to the second
        if($signDx[0] <= 0) {
            // Want alternating signs
            $ii = -1;
            if($signDx[0] == $signDx[1]) {
                for($z=1; $z<$lenX; $z++) {
                    $x[$z-1] = $x[$z];
                    $index[$z-1] = $index[$z];
                }
                $lenX--;
                // Delete last value
                unset($x[$lenX]);
                unset($index[$lenX]);
            }
        }
        else {
            // First point is smaller than the second
            $ii = 0;
            if($signDx[0] == $signDx[1]) {
                array_pop($x);
                array_pop($index);
                $lenX--;
            }
        }
        // Preallocate max number of maxima
        $maxPeaks = ceil($lenX/2);
        $peakLoc = zeros($maxPeaks);
        $peakMag = zeros($maxPeaks);
        $cInd = 0;
        $foundPeak = false;
        $tempLoc = 0;
        // Loop through extrema which should be peaks and then valleys
        while(($ii+1) < $lenX) {
            $ii++;
            // This is a peak
            // Reset peak finding if we had a peak and the next peak is bigger
            // than the last or the left min was small enough to reset
            if($foundPeak) {
                $tempMag = $minMag;
                $foundPeak = false;
            }
            // Make sure we don't iterate past the length of our vector
            if($ii == ($lenX-1)) {
                // We assign the last point differently out of the loop
                break;
            }
            // Found new peak that was lager than temp mag and selectivity larger than the minimum to its left
            if(($x[$ii] > $tempMag) && ($x[$ii] > ($leftMin + $sel))) {
                $tempLoc = $ii;
                $tempMag = $x[$ii];
            }
            //$ii++;
            // Move onto the valley
            // Come down at least sel from peak
            //print((($foundPeak==true)?"true":"false") ." ".$tempMag." ".$x[$ii]."<br>");
            if((!$foundPeak) && ($tempMag > ($sel + $x[$ii]))) {
                $foundPeak = true; // We have found a peak
                $leftMin = $x[$ii];
                $peakLoc[$cInd] = $tempLoc; // Add peak to index
                $peakMag[$cInd] = $tempMag;
                $cInd++;
            }
            else if($x[$ii] < $leftMin) {
                // New left minima
                $leftMin = $x[$ii];
            }
        }
        // Check end point
        if(($x[$lenX-1] > $tempMag) && ($x[$lenX-1] > $leftMin + $sel)) {
            $peakLoc[$cInd] = $lenX-1;
            $peakMag[$cInd] = $x[$lenX-1];
            $cInd++;
        }
        else if(!$foundPeak && ($tempMag > $minMag)) {
            // Check if we still need to add the last point
                $peakLoc[$cInd] = $tempLoc;
                $peakMag[$cInd] = $tempMag;
                $cInd++;
        }
        // Create output
        $peakInds = array();
        $peakMags = array();
        for($i=0; $i<$cInd; $i++) {
            $peakInds[] = $index[$peakLoc[$i]];
            $peakMags[] = $peakMag[$i];
        }
    }
    else {
        // This is a monotone function where an endpoint is the only peak
        $maxIndex = maxIndex($x);
        $peakMags = $maxIndex['value'];
        $xInd = $maxIndex['index'];
        if($peakMags > $minMag + $sel) {
            $peakInds[] = $index[$xInd];
        }
        else {
            unset($peakMags);
            unset($peakLoc);
        }
    }
    // Apply threshold value.  Since always finding maxima it will always be larger than the thresh.
    if($thresh!=null) {
        $peakIndsNew = array();
        $peakMagsNew = array();
        for($i=0; $i<count($thresh); $i++) {
            if($peakMags[$i]>$thresh) {
                $peakIndsNew[] = $peakInds[$i];
                $peakMagsNew[] = $peakMags[$i];
            }
        }
        $peakInds = $peakIndsNew;
        $peakMags = $peakMagsNew;
    }
    // Change sign of data if was finding minima
    if($extrema < 0) {
        $peakMags = mulArray($peakMags, -1);
        $x0 = mulArray($x0, -1);
    }
    $response = array();
    $response['values'] = $peakMags;
    $response['indexs'] = $peakInds;
    return $response;
}

?>