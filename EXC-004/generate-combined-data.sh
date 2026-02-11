#!/bin/bash

## This script processes metagenomic sequencing data from multiple cultures
## It organizes genome bins (FASTA files) and metadata into a combined directory

## ==============================================================================
## STEP 1: Create the COMBINED-DATA directory
## ==============================================================================
## This creates a new folder at the same level as RAW-DATA to store all results
echo "Step 1: Creating COMBINED-DATA directory..."
mkdir -p COMBINED-DATA
echo "✓ COMBINED-DATA directory created"
echo ""

## ==============================================================================
## STEP 2: Load the sample translation file from RAW-DATA directory
## ==============================================================================
## This file maps library names (DNA57, DNA58, etc.) to culture names (CO64, CO83, etc.)
## The file is located inside the RAW-DATA directory
echo "Step 2: Loading sample translation mappings..."

## Check if the RAW-DATA directory exists
if [ ! -d "RAW-DATA" ]; then
    echo "ERROR: RAW-DATA directory not found!"
    echo "Please make sure you are running this script from the correct directory."
    exit 1
fi

## Check if the translation file exists inside RAW-DATA
if [ ! -f "RAW-DATA/sample-translation.txt" ]; then
    echo "ERROR: sample-translation.txt not found in RAW-DATA directory!"
    exit 1
fi

## Create an associative array to store library->culture mappings
declare -A CULTURE_MAP

## Read the file and populate the mapping
## We skip the header line and process each data line
while IFS=$'\t' read -r library culture rest; do
    ## Skip empty lines and header
    if [ -n "$library" ] && [ "$library" != "library" ]; then
        CULTURE_MAP[$library]=$culture
        echo "  Mapped: $library -> $culture"
    fi
done < RAW-DATA/sample-translation.txt

echo "✓ Loaded ${#CULTURE_MAP[@]} sample mappings"
echo ""

## ==============================================================================
## STEP 3: Process each directory in RAW-DATA
## ==============================================================================
echo "Step 3: Processing each sample directory..."
echo ""

## Loop through all directories in RAW-DATA
for dir in RAW-DATA/*/; do
    ## Remove trailing slash and get just the directory name (e.g., DNA57)
    library=$(basename "$dir")
    
    echo "----------------------------------------"
    echo "Processing library: $library"
    echo "----------------------------------------"
    
    ## Get the culture name from our mapping
    culture=${CULTURE_MAP[$library]}
    
    ## Skip if we don't have a mapping for this library
    if [ -z "$culture" ]; then
        echo "WARNING: No culture mapping found for $library, skipping..."
        echo ""
        continue
    fi
    
    echo "Culture name: $culture"
    
    ## ----------------------------------------------------------------------
    ## STEP 3a: Copy checkm.txt and gtdb.gtdbtk.tax files
    ## ----------------------------------------------------------------------
    echo ""
    echo "Copying metadata files..."
    
    ## Copy CheckM results (genome completeness/contamination info)
    ## These files are in the library directory (e.g., RAW-DATA/DNA57/)
    if [ -f "${dir}checkm.txt" ]; then
        cp "${dir}checkm.txt" "COMBINED-DATA/${culture}-CHECKM.txt"
        echo "  ✓ Copied checkm.txt -> ${culture}-CHECKM.txt"
    else
        echo "  WARNING: checkm.txt not found for $library"
    fi
    
    ## Copy GTDB taxonomy results
    if [ -f "${dir}gtdb.gtdbtk.tax" ]; then
        cp "${dir}gtdb.gtdbtk.tax" "COMBINED-DATA/${culture}-GTDB-TAX.txt"
        echo "  ✓ Copied gtdb.gtdbtk.tax -> ${culture}-GTDB-TAX.txt"
    else
        echo "  WARNING: gtdb.gtdbtk.tax not found for $library"
    fi
    
    ## ----------------------------------------------------------------------
    ## STEP 3b: Load CheckM data to determine MAG vs BIN classification
    ## ----------------------------------------------------------------------
    echo ""
    echo "Loading CheckM data for quality assessment..."
    
    ## Create associative arrays to store completion and contamination values
    declare -A COMPLETION
    declare -A CONTAMINATION
    
    ## Read the checkm.txt file from the library directory
    if [ -f "${dir}checkm.txt" ]; then
        ## Skip the header lines and read bin quality data
        while IFS=$'\t' read -r bin_id marker_lineage num_genomes num_markers num_marker_sets comp contam strain_het; do
            ## Skip header lines (they contain text like "Bin Id")
            if [[ "$bin_id" != "Bin Id"* ]] && [ -n "$bin_id" ]; then
                COMPLETION[$bin_id]=$comp
                CONTAMINATION[$bin_id]=$contam
                echo "  $bin_id: Completion=${comp}%, Contamination=${contam}%"
            fi
        done < "${dir}checkm.txt"
    fi
    
    ## ----------------------------------------------------------------------
    ## STEP 3c: Process FASTA files in the bins/ directory
    ## ----------------------------------------------------------------------
    echo ""
    echo "Processing FASTA files..."
    
    ## Initialize counters for sequential numbering
    ## MAGs and BINs are numbered independently
    mag_counter=1
    bin_counter=1
    
    ## Check if bins directory exists inside the library directory
    if [ ! -d "${dir}bins" ]; then
        echo "  WARNING: bins/ directory not found for $library"
        echo ""
        continue
    fi
    
    ## Process each FASTA file in the bins directory (e.g., RAW-DATA/DNA57/bins/)
    for fasta in "${dir}bins/"*.fasta; do
        ## Skip if no files match the pattern
        [ -e "$fasta" ] || continue
        
        ## Get the base filename without path or extension
        bin_name=$(basename "$fasta" .fasta)
        
        echo ""
        echo "  Processing: $bin_name"
        
        ## --------------------------------------------------------------
        ## SPECIAL CASE: Handle unbinned sequences
        ## --------------------------------------------------------------
        if [[ "$bin_name" == *"unbinned"* ]]; then
            ## Output file goes in COMBINED-DATA directory
            output_file="COMBINED-DATA/${culture}_UNBINNED.fa"
            echo "    Type: UNBINNED sequences"
            echo "    Output: $output_file"
            
            ## Copy and modify FASTA headers to include culture name AND bin identifier
            ## This ensures all sequences have unique identifiers across all files
            ## Format: >CULTURE_BINNAME_originalheader
            awk -v culture="$culture" -v bin="UNBINNED" '{
                if ($0 ~ /^>/) {
                    ## This is a header line (starts with >)
                    ## Add culture and bin prefix to make it unique
                    print ">" culture "_" bin "_" substr($0, 2)
                } else {
                    ## This is a sequence line, print as-is
                    print $0
                }
            }' "$fasta" > "$output_file"
            
            echo "    ✓ Created with unique deflines"
            continue
        fi
        
        ## --------------------------------------------------------------
        ## REGULAR BINS/MAGs: Determine classification
        ## --------------------------------------------------------------
        
        ## Get quality metrics for this bin from CheckM data
        comp=${COMPLETION[$bin_name]}
        contam=${CONTAMINATION[$bin_name]}
        
        ## Determine if this is a MAG or BIN
        ## MAG criteria: completion ≥ 50% AND contamination < 5%
        is_mag=0
        if [ -n "$comp" ] && [ -n "$contam" ]; then
            ## Use bc for floating point comparison
            if (( $(echo "$comp >= 50" | bc -l) )) && (( $(echo "$contam < 5" | bc -l) )); then
                is_mag=1
            fi
        fi
        
        ## Set output filename based on classification
        ## All output files go to COMBINED-DATA directory
        if [ $is_mag -eq 1 ]; then
            ## This is a MAG (Metagenome-Assembled Genome) - high quality
            seq_num=$(printf "%03d" $mag_counter)
            output_file="COMBINED-DATA/${culture}_MAG_${seq_num}.fa"
            echo "    Type: MAG (High quality)"
            echo "    Quality: Completion=${comp}%, Contamination=${contam}%"
            echo "    Number: MAG_$seq_num"
            mag_counter=$((mag_counter + 1))
        else
            ## This is a BIN (lower quality)
            seq_num=$(printf "%03d" $bin_counter)
            output_file="COMBINED-DATA/${culture}_BIN_${seq_num}.fa"
            echo "    Type: BIN (Standard quality)"
            if [ -n "$comp" ]; then
                echo "    Quality: Completion=${comp}%, Contamination=${contam}%"
            else
                echo "    Quality: No CheckM data available"
            fi
            echo "    Number: BIN_$seq_num"
            bin_counter=$((bin_counter + 1))
        fi
        
        echo "    Output: $output_file"
        
        ## --------------------------------------------------------------
        ## Modify FASTA headers to include culture name AND bin type/number for uniqueness
        ## Read from RAW-DATA, write to COMBINED-DATA
        ## Format: >CULTURE_BINTYPE_NUMBER_originalheader (e.g., >CO64_MAG_001_contig_1)
        ## --------------------------------------------------------------
        if [ $is_mag -eq 1 ]; then
            bin_identifier="MAG_${seq_num}"
        else
            bin_identifier="BIN_${seq_num}"
        fi
        
        awk -v culture="$culture" -v bin_id="$bin_identifier" '{
            if ($0 ~ /^>/) {
                ## This is a header line (starts with >)
                ## Add culture and bin identifier prefix to make it unique across all samples and bins
                print ">" culture "_" bin_id "_" substr($0, 2)
            } else {
                ## This is a sequence line, print as-is
                print $0
            }
        }' "$fasta" > "$output_file"
        
        echo "    ✓ Created with unique deflines"
    done
    
    ## Clean up the associative arrays for the next iteration
    unset COMPLETION
    unset CONTAMINATION
    
    echo ""
    echo "✓ Completed processing $library ($culture)"
    echo ""
done

## ==============================================================================
## STEP 4: Summary
## ==============================================================================
echo "========================================"
echo "PROCESSING COMPLETE!"
echo "========================================"
echo ""
echo "Summary of results:"
echo "  Total files in COMBINED-DATA: $(ls -1 COMBINED-DATA | wc -l)"
echo ""
echo "Files created:"
ls -1 COMBINED-DATA | head -20
if [ $(ls -1 COMBINED-DATA | wc -l) -gt 20 ]; then
    echo "  ... and $(( $(ls -1 COMBINED-DATA | wc -l) - 20 )) more files"
fi
echo ""
echo "✓ All done! Check the COMBINED-DATA directory for results."