
for file in test_1000_samples/*
do
    nlines=$(wc -l $file|awk '{print $1}')
    #echo $nlines
    if [ "$nlines" -eq 1 ]; then
    echo $file
    fi
done
