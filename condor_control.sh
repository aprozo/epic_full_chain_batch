#!/bin/bash

#change the following parameters
runningTimeLimitHours=48 #hours
sleepTime=10             #seconds
#  ========================================================
username=$(whoami)                                 # Set the username for which to check the jobs
runningTime=0                                      # Initialize the running time to 0
condorQFile='condor.out'                           # File to save the output of condor_q
runningTimeLimit=$((runningTimeLimitHours * 3600)) #seconds
flag=1                                             # Flag to check if it is the first time

# Function to check if all jobs for the user are finished
check_jobs() {
    condor_q >$condorQFile              # Save the output of condor_q to a file
    head_jobs=$(head -n 5 $condorQFile) # Get first 5 lines
    # Extract the line containing the user's job information
    userJobsLine=$(grep "Total for $username" $condorQFile)
    # Extract the number of jobs from the user's line
    leftJobs=$(echo "$userJobsLine" | awk '{print $4}')     # 4th column (word) in the line
    idleJobs=$(echo "$userJobsLine" | awk '{print $10}')    # 10th column (word) in the line
    runningJobs=$(echo "$userJobsLine" | awk '{print $12}') # 12th column (word) in the line
    heldJobs=$(echo "$userJobsLine" | awk '{print $14}')    # 14th column (word) in the line

    #example userJobLine:
    #Total for username: 10 jobs; 1 completed, 2 removed, 3 idle, 0 running, 4 held, 0 suspended
    if [ "$flag" = "1" ]; then
        totalJobs=$leftJobs
        flag=0
    fi

    # Check if all jobs are finished
    if [ "$leftJobs" = "0" ]; then
        echo "Jobs for user $username are finished! -------------------> DONE"
        return 1
    # Check if the running time exceeds the limit
    elif [ $runningTime -gt $runningTimeLimit ]; then
        echo "Jobs for user $username are running for too long! Removing all remaining jobs..."
        condor_rm $username
        echo "Further resubmission needed"
        return 1
    # Check if there are held jobs
    elif [ "$heldJobs" != "0" ] && [ "$idleJobs" = "0" ] && [ "$runningJobs" = "0" ]; then
        echo "There are held jobs for user $username. Removing all remaining jobs..."
        condor_rm $username
        echo "Further resubmission needed"
        return 1

    else
        # If there are still running jobs, wait for 10 seconds
        # Increment runningTime by sleep interval
        runningTime=$((runningTime + sleepTime))
        # Calculate percentages and time metrics
        percentage_time=$((runningTime * 100 / runningTimeLimit))
        percentage_done=$((100 - $leftJobs * 100 / totalJobs))
        runningTimeHours=$((runningTime / 3600))

        # Print status on a single line with carriage returns to overwrite previous output
        printf "Jobs: %d/%d (%d%% done) | Time: %dh %dm %ds (%d%% of limit) | Waiting %ds...\r" \
            $((totalJobs - $leftJobs)) $totalJobs $percentage_done \
            $((runningTime / 3600)) $(((runningTime % 3600) / 60)) $((runningTime % 60)) \
            $percentage_time $sleepTime

        # Sleep for the specified interval
        sleep $sleepTime
    fi
}
# Main loop to continuously check job status
while :; do
    check_jobs || break # If all jobs are finished, exit the loop
done

rm $condorQFile
