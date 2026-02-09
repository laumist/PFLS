#skript to get the output for excercise 3

for file in $@   #grabs all arguments
do
    #define variables
    no_seq=$(grep '>' $file | wc -l | awk '{print $1}')
    length=$(awk '!/>/{total += gsub(/[AaTtGgCc]/, "")} END{print total}' $file) #counts A,a,T,t,G,g;C,c over the whole document, excludes lines starting with '>'
    #counts all elements of sequence
    seq_length=$(awk '
                    />/{
                        if (len) print len; len=0; next
                    }
                    !/>/{len = len + length}
                    END{print len}
                    ' $file ) 
    
    CG=$(awk '!/>/ {gc_count += gsub(/[GgCc]/, "")} END {print gc_count}' $file)
    line=$(awk '!/>/{print $0}' $file)

    #just calculate if the sequence has just letters
    if ! grep -v '^>' "$file" | grep -q '[^A-Za-z ]'; then
        #echo what meren wants to see
        echo "FASTA File Statistics:" #header 
        echo "----------------------"
        echo "Number of sequences: $no_seq" 
        echo "Total length of sequences: $length"
        echo "Length of the longest sequence:" $(echo "$seq_length" | sort -n | tail -n 1)
        echo "Length of the shortest sequence:" $(echo "$seq_length" | sort -nr | tail -n 1)
        echo "Average sequence length:" $(echo "scale=3; $length / $no_seq" | bc)
        echo "GC Content (%):" $(echo "scale=3; $CG * 100 / $length" | bc)
  
    else 
        echo "FASTA File Statistics:"
        echo "----------------------"
        echo "Number of sequences: 0"
        echo "Total length of sequences: 0"
        echo "Length of the longest sequence: 0"
        echo "Length of the shortest sequence: 0"
        echo "Average sequence length: 0"
        echo "GC Content (%): 0"
    fi 
done
