#!/bin/bash

# Absolute paths (no alias)
LGEN="/usr/local/bin/latticegen -randseed time"
FPL="/usr/local/bin/fplll -d 0.75"
LLL="./lll_exec"

BITS=50
DIMS=(5 10 50 100)
TRIALS=10

for DIM in "${DIMS[@]}"; do
    echo "==== Testing dimension ${DIM} x ${DIM} ===="

    for ((i = 1; i <= TRIALS; i++)); do
        echo "→ Trial $i"

        # Step 1: Write dimension header to input
        echo "$DIM $DIM" > clean

        # Step 2: Generate random lattice and clean brackets
        $LGEN u $DIM $BITS > file
        sed 's/[][]//g' file >> clean

        # Step 3: Run your implementation and measure time
        echo -n "   Your LLL time:     "
        /usr/bin/time -p $LLL < clean > output 2> time_lll.log
        grep user time_lll.log | awk '{print $2 " sec"}'

	# Step 4: Run fplll and measure time
        echo -n "   FPLLL time:        "
        /usr/bin/time -p $FPL file > outfpl 2> time_fpl.log   
        grep user time_fpl.log | awk '{print $2 " sec"}'



        # Step 5: Clean fplll output
        sed 's/[][]//g' outfpl > cleanfpl

        # Step 6: Count differing vectors
        diff_count=$(diff output cleanfpl | grep '^<' | wc -l)
        echo "   Vectors differing from FPLLL:  $diff_count"

        # Step 7: Run internal verifier
        echo -n "   Verifier:          "
        if $LLL < clean | grep -q "Failure"; then
            echo "Fail"
        else
            echo "Pass"
        fi

        echo ""
    done
done

