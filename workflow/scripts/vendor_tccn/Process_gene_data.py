import numpy as np

def process_gene_data(filename):
    """
    Process gene data from a file and return arrays and gene_dict.

    Parameters:
    filename (str): The path to the input file.

    Returns:
    tuple: A tuple containing arrays and gene_dict.
    """
    # Read the file
    with open(filename, 'r') as file:
        lines = file.readlines()

    if len(lines) < 4:
        return 'The file does not contain enough lines.'

    # Find the values in the second column of the third and fourth lines
    third_line_second_column = int(lines[2].split()[1])
    fourth_line_second_column = int(lines[3].split()[1])

    # If the order is descending, reorder the exons for each transcript
    if third_line_second_column > fourth_line_second_column:
        
        # Initialize variables
        new_lines = []  # Store the processed content
        current_transcript_exons = []  # Store the exon lines of the current transcript
        transcript_found = False  # Flag to indicate if a transcript is found

        # Iterate through the file content
        for line in lines:
            # Check if the current line is a transcript
            if line.startswith('transcript'):
                if transcript_found:
                    
                    # If a transcript was found before, reverse the order of exons and add to new_lines
                    new_lines.extend(current_transcript_exons[::-1])
                    
                # Reset variables
                current_transcript_exons = []
                transcript_found = True
                new_lines.append(line)  # Directly add the transcript line to new_lines

            # Check if the current line is an exon
            elif 'exon' in line and transcript_found:
                current_transcript_exons.append(line)
                
            # Add other lines directly to new_lines
            else:
                new_lines.append(line)

        # Add a newline character to the last element
        current_transcript_exons[-1] += '\n'

        # Add the exons of the last transcript
        new_lines.extend(current_transcript_exons[::-1])

        # Write the processed content to a new file
        with open(filename, 'w') as file:
            file.writelines(new_lines)

    # Initialize the array
    exon_starts = []

    # Read the file and process
    with open(filename, 'r') as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        
        # Check if the line is an exon line
        if line.startswith('exon'):
            
            # Get the value in the third column
            third_column = int(line.split()[2])
            
            # Check if the next line is a transcript
            if i < len(lines) - 1 and lines[i + 1].startswith('transcript') or i == len(lines) - 1:
                
                # If the next line is a transcript, get the value in the second column
                second_column = int(line.split()[1])
                exon_starts.append(second_column)
            else:
                exon_starts.append(third_column)

    exon_unique = sorted(set(exon_starts))

    # Assume exon_starts_unique_sorted is the processed list
    exon_starts_unique_sorted = exon_unique

    # Initialize the dictionary
    exon_lengths = {i + 1: None for i in range(len(exon_starts_unique_sorted))}

    with open(filename, 'r') as file:
        lines = file.readlines()

    for i, start in enumerate(exon_starts_unique_sorted):
        max_length = None  # Store the maximum length difference for the current number
        for line in lines:
            if 'exon' in line:
                parts = line.split()
                second_column = int(parts[1])
                third_column = int(parts[2])
                if second_column == start or third_column == start:
                    length_difference = third_column - second_column
                    # Update the maximum length difference
                    if max_length is None or length_difference > max_length:
                        max_length = length_difference
                        
        # Store the maximum length difference in the dictionary
        exon_lengths[i + 1] = max_length

    gene_dict = exon_lengths

    transcript_exon_positions = []  # Store the exon positions for each transcript
    current_transcript_positions = []  # Store the exon positions for the current transcript

    with open(filename, 'r') as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if line.startswith('transcript'):
            
            # If not the first transcript, save the exon positions of the previous transcript
            if current_transcript_positions:
                transcript_exon_positions.append(current_transcript_positions)
            current_transcript_positions = []  # Reset the exon positions list for the current transcript
        elif line.startswith('exon'):
            parts = line.split()
            second_column = int(parts[1])
            third_column = int(parts[2])
            
            # Check if the next line is a transcript or the file has ended
            if i < len(lines) - 1 and lines[i + 1].startswith('transcript') or i == len(lines) - 1:
                # For the last exon, check the second element
                position = exon_starts_unique_sorted.index(second_column) + 1
            else: 
                # For regular exons, check the third element
                position = exon_starts_unique_sorted.index(third_column) + 1
            current_transcript_positions.append(position)

    # Add the exon positions of the last transcript
    if current_transcript_positions:
        transcript_exon_positions.append(current_transcript_positions)

    arrays = transcript_exon_positions
    return arrays, gene_dict