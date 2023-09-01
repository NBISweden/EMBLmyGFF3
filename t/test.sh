#!/usr/bin/env bash

## SETUP VIRTUAL ENVIRONMENT FOR TESTING
if [ "x$TRAVIS" != "x" ]; then
    if [ -d t/venv ]; then
        rm -rf t/venv
    fi

    virtualenv t/venv
    source t/venv/bin/activate
    trap 'deactivate' EXIT
fi

python setup.py install


## RUN TESTS
cd t

# Remove line related to dates
if [ "$(uname)" == "Darwin" ]; then
    sed -Ei '' -e '/^DT   .*/d' \
        -e '/^RL   .*/d' \
       *.embl  
else
    sed -Ei -e '/^DT   .*/d' \
        -e '/^RL   .*/d' \
       *.embl 
fi

SUCCESS=0
FAIL=0

for NAME in augustus maker prokka prokka_disorder dbxref_test aa; do
    RESULT_FILE="EMBLmyGFF3-${NAME}-example.embl"
    EXPECTED_FILE="EMBLmyGFF3-${NAME}-test.embl"
    cp $EXPECTED_FILE $EXPECTED_FILE.copy
    [ -f "$RESULT_FILE" ] && rm $RESULT_FILE
    ../examples/${NAME}_example.py

    # Remove line related to dates
    if [ "$(uname)" == "Darwin" ]; then
        sed -Ei '' -e '/^DT   .*/d' \
            -e '/^RL   .*/d' \
        $RESULT_FILE $EXPECTED_FILE.copy
    else
        sed -Ei -e '/^DT   .*/d' \
            -e '/^RL   .*/d' \
        $RESULT_FILE $EXPECTED_FILE.copy
    fi

    if diff -q "$RESULT_FILE" "$EXPECTED_FILE.copy"; then
        SUCCESS=$(( $SUCCESS + 1 ))
    else
        diff "$RESULT_FILE" "$EXPECTED_FILE.copy"
        FAIL=$(( $FAIL + 1 ))
    fi
    
    [ -f "$RESULT_FILE" ] && rm $RESULT_FILE
    [ -f "$EXPECTED_FILE.copy" ] && rm $EXPECTED_FILE.copy

done

if [ $FAIL -eq 0 ]; then
    echo "All tests successfull"
    exit 0
fi

echo "Failed $FAIL out of 6 tests"
exit 1
