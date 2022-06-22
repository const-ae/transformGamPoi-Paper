
# Functions from https://unix.stackexchange.com/a/216475
# initialize a semaphore with a given number of tokens
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

# run the given command asynchronously and pop/push tokens
run_with_lock(){
    local x
    # this read waits until there is something to read
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    # push the return code of the command to the semaphore
    printf '%.3d' $? >&3
    )&
}

N=1
open_sem $N
for notebook in notebooks/*.Rmd; do
  run_with_lock R -e "rmarkdown::render('$notebook', output_format = 'html_document')"
done

wait
echo "Done rendering"

rm notebooks/*-tikzDictionary

echo "Done"