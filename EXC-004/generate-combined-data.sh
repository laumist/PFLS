#!/bin/bash

echo "Step 1: Creating COMBINED-DATA directory..."
mkdir -p COMBINED-DATA
echo "✓ COMBINED-DATA directory created"
echo ""

echo "Step 2: Loading sample translation mappings..."

if [ ! -d "RAW-DATA" ]; then
    echo "ERROR: RAW-DATA directory not found!"
    exit 1
fi

if [ ! -f "RAW-DATA/sample-translation.txt" ]; then
    echo "ERROR: sample-translation.txt not found in RAW-DATA directory!"
    exit 1
fi

declare -A CULTURE_MAP

while IFS=$'\t' read -r library culture rest; do
    if [ -n "$library" ] && [ "$library" != "library" ]; then
        CULTURE_MAP[$library]=$culture
        echo "  Mapped: $library -> $culture"
    fi
done < RAW-DATA/sample-translation.txt

echo "✓ Loaded ${#CULTURE_MAP[@]} sample mappings"
echo ""

echo "Step 3: Processing each sample directory..."
echo ""

for dir in RAW-DATA/*/; do

    library=$(basename "$dir")
    echo "----------------------------------------"
    echo "Processing library: $library"
    echo "----------------------------------------"

    culture=${CULTURE_MAP[$library]}

    if [ -z "$culture" ]; then
        echo "WARNING: No culture mapping found for $library, skipping..."
        echo ""
        continue
    fi

    echo "Culture name: $culture"
    echo ""
    echo "Copying metadata files..."

    if [ -f "${dir}checkm.txt" ]; then
        cp "${dir}checkm.txt" "COMBINED-DATA/${culture}-CHECKM.txt"
        echo "  ✓ Copied checkm.txt -> ${culture}-CHECKM.txt"
    else
        echo "  WARNING: checkm.txt not found for $library"
    fi

    if [ -f "${dir}gtdb.gtdbtk.tax" ]; then
        cp "${dir}gtdb.gtdbtk.tax" "COMBINED-DATA/${culture}-GTDB-TAX.txt"
        echo "  ✓ Copied gtdb.gtdbtk.tax -> ${culture}-GTDB-TAX.txt"
    else
        echo "  WARNING: gtdb.gtdbtk.tax not found for $library"
    fi

    echo ""
    echo "Processing FASTA files..."

    mag_counter=1
    bin_counter=1

    if [ ! -d "${dir}bins" ]; then
        echo "  WARNING: bins/ directory not found for $library"
        echo ""
        continue
    fi

    for fasta in "${dir}bins/"*.fasta; do
        [ -e "$fasta" ] || continue

        bin_name=$(basename "$fasta" .fasta)

        echo ""
        echo "  Processing: $bin_name"

        if [[ "$bin_name" == *"unbinned"* ]]; then
            output_file="COMBINED-DATA/${culture}_UNBINNED.fa"

            awk -v culture="$culture" -v bin="UNBINNED" '{
                if ($0 ~ /^>/) {
                    print ">" culture "_" bin "_" substr($0, 2)
                } else {
                    print $0
                }
            }' "$fasta" > "$output_file"

            continue
        fi

        comp=""
        contam=""

        # Use awk to read completeness ($13) and contamination ($14) from checkm.txt
        if [ -f "${dir}checkm.txt" ]; then
            read comp contam < <(awk -v bin="$bin_name" '$1 ~ bin && $1 != "Bin" {print $13, $14; exit}' "${dir}checkm.txt")
        fi

        is_mag=0
        if [ -n "$comp" ] && [ -n "$contam" ]; then
            if (( $(echo "$comp >= 50" | bc -l) )) && (( $(echo "$contam < 5" | bc -l) )); then
                is_mag=1
            fi
        fi

        if [ $is_mag -eq 1 ]; then
            seq_num=$(printf "%03d" $mag_counter)
            output_file="COMBINED-DATA/${culture}_MAG_${seq_num}.fa"
            mag_counter=$((mag_counter + 1))
        else
            seq_num=$(printf "%03d" $bin_counter)
            output_file="COMBINED-DATA/${culture}_BIN_${seq_num}.fa"
            bin_counter=$((bin_counter + 1))
        fi

        if [ $is_mag -eq 1 ]; then
            bin_identifier="MAG_${seq_num}"
        else
            bin_identifier="BIN_${seq_num}"
        fi

        awk -v culture="$culture" -v bin_id="$bin_identifier" '{
            if ($0 ~ /^>/) {
                print ">" culture "_" bin_id "_" substr($0, 2)
            } else {
                print $0
            }
        }' "$fasta" > "$output_file"

    done

    echo ""
    echo "✓ Completed processing $library ($culture)"
    echo ""
done

echo "========================================"
echo "PROCESSING COMPLETE!"
echo "========================================"