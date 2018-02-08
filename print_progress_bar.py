def print_progress_bar(current_status, completion_status, message):
    """
    Prints the progress bar for the current task given its current progress,
    the total completion state (all abstracted to a number), and a message to
    display next to the bar.

    USAGE:
        current_status and completion_status should be numbers.
        message should be a string.
    """
    import sys
    bar_length = 50
    percent = float((current_status) / completion_status)
    hashes = '#' * int(round(percent * bar_length))
    spaces = ' ' * (bar_length - len(hashes))
    sys.stdout.write("\r" + message + ": [{0}] {2} \t {1}%".format(hashes
                     + spaces, int(round(percent * 100)), str(current_status)
                     + '/' + str(completion_status)))
    sys.stdout.flush()
