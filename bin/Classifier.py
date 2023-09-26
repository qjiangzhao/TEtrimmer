import os
import sys
import repeatmodeler_classify  # Import the repeatmodeler_classify module

if __name__ == "__main__":
    # Check if the user provided the path to the consensus file as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python Classifier.py <consensus_file>")
        sys.exit(1)

    # Get the consensus file path from the command-line argument
    consensus_file = sys.argv[1]

    # Call repeatmodeler_classify.classify with the provided consensus file path
#   "RepeatClassifier", "-consensi", consensus_file
    repeatmodeler_classify.classify(consensus_file)