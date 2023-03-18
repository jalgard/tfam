'''
    Terminal .FAsta Mutator

    written by Evgeny Gerasimov, 2023
    Lomonosov Moscow State University, Biological Faculty, 330

'''

import sys, random, os, argparse, time
from collections import defaultdict

# try importing termcolor module
has_termcolor = True
try:
    from termcolor import colored, cprint

except ImportError:
    has_termcolor = False

font_styles = {
    'text_fg'         : 'black',
    'text_bg'         : None,
    'text_attrs'      : [],
    'info_fg'         : 'green',
    'info_bg'         : None,
    'info_attrs'      : [],
    'highlight_fg'    : 'white',
    'highlight_bg'    : 'on_red',
    'highlight_attrs' : ['bold'],
    'marked_fg'       : 'red',
    'marked_bg'       : None,
    'marked_attrs'    : ['bold'],
    'motif_fg'        : 'yellow',
    'motif_bg'        : 'on_blue',
    'motif_attrs'     : ['bold']
}

default_output_dest = sys.stdout

'''
    WRITERS
'''
def printer(what, message_type = 'text', newline = True, stream = default_output_dest, has_termcolor = has_termcolor, styles = font_styles):
    """ general function to print output """

    if has_termcolor:
        what = colored(what, color = styles['{}_fg'.format(message_type)], on_color = styles['{}_bg'.format(message_type)], attrs = styles['{}_attrs'.format(message_type)])

    line_ending = '\n'
    if newline == False:
        line_ending = ''

#    elif newline == 'Windows':
#        line_ending = '\r\n'

    stream.writelines('{}{}'.format(what, line_ending))


def write_fasta(fasta_entries, stream = default_output_dest):
    """ default output of processed data """

    for i in range(len(fasta_entries)):
        printer('>', message_type = 'highlight', newline = False, stream=stream)
        printer(fasta_entries[i][0], message_type = 'info', stream=stream)
        printer(fasta_entries[i][1], message_type = 'text', stream=stream)


'''
    READERS
'''
def fasta_reader(fasta_handle):
    """ reads .fasta into list """

#    Fast version of .fasta reader uses Tanya's .replace() hack to remove \r, \n
#    from the sequence string and join() to create initial sequence string
#    from list
#    Function reads fasta entries as is, does not split anything

    fasta_lines  = fasta_handle.readlines()
    # add dummy entry at the end of the data to trigger last round of
    # sequence_accu insertion
    fasta_lines += ['>\n']
    fasta_entries = []
    sequence_accu = []
    for line in fasta_lines:
        if line[0] == '>':
            if len(sequence_accu) > 0:
                # fast sequence string formatting
                fasta_entries[-1][1] = ''.join(sequence_accu).replace('\n', '').replace('\r', '')
                sequence_accu = []

            fasta_entries.append([line[1:].strip(), ''])

        else:
            sequence_accu.append(line)

    # return everything except dummy entry
    return fasta_entries[:-1]


def load_filter_dict(input_file):
    """ reads list file for filtering """

    filter_dict = {}
    with open(input_file, 'r') as ifile:
        for line in ifile:
            filter_item = line.strip()
            if filter_item[0] == '>':
                filter_item = filter_item[1:]

            filter_dict[filter_item] = 1

    return filter_dict


'''
    SUPPORTED OPERATIONS SECTION
'''
def EntryUppercase(input_entries, **options):
    """ all sequence uppercase """

    for i in range(len(input_entries)):
        input_entries[i][1] = input_entries[i][1].upper()

def EntryLowercase(input_entries, **options):
    """ all sequence lowercase """

    for i in range(len(input_entries)):
        input_entries[i][1] = input_entries[i][1].lower()

def ListRemove(input_entries, **options):
    """ remove listed sequences """

    filter = options.get('filter_dict')
    filtered_entries = []
    for i in range(len(input_entries)):
        if input_entries[i][0] not in filter:
            filtered_entries.append(input_entries[i])

    input_entries[:] = list(filtered_entries)

def ListKeep(input_entries, **options):
    """ keep listed sequences """

    filter = options.get('filter_dict')
    filtered_entries = []
    for i in range(len(input_entries)):
        if input_entries[i][0] in filter:
            filtered_entries.append(input_entries[i])

    input_entries[:] = list(filtered_entries)

def InfoN50(input_entries, **options):
    """ only output N50 of .fasta provided """

    total_len = 0
    for i in range(len(input_entries)):
        input_entries[i].append(len(input_entries[i][1]))
        total_len += input_entries[i][-1]

    running_sum = 0
    for i in sorted(input_entries, key = lambda x : x[2], reverse = True):
        running_sum += i[2]
        if running_sum > total_len / 2:
            printer('N50\t', message_type = 'highlight', newline = False)
            printer(i[2], newline = False)
            printer('\t', newline = False)
            printer(i[0])
            # EXIT point, nothing to be done
            exit(0)


def InfoEntryLength(input_entries, **options):
    """ print length of each entry as table """

    for i in range(len(input_entries)):
        printer(input_entries[i][0], newline = False, message_type = 'text')
        printer('\t', newline = False, message_type = 'text')
        printer(len(input_entries[i][1]), newline = True, message_type = 'highlight')
    # EXIT point, nothing to be done
    exit(0)

ACTIONS = {
'upper'  :    EntryUppercase,          # action:  upper
'lower'  :    EntryLowercase,          # action:  lower
'remove' :    ListRemove,              # action:  remove     [needs --list]
'keep'   :    ListKeep,                # action:  keep       [needs --list]
#'rename' :    EntryRename,
'N50'    :    InfoN50,
'len'    :    InfoEntryLength,
#'rc'     :    EntryRecvom
}





'''
    ARGUMENTS PARSER
'''
TfamOptionsParser = argparse.ArgumentParser()
TfamOptionsParser.add_argument('--in', dest='fasta', action='store', help="Input .fasta file name (or stdin if not set")
TfamOptionsParser.add_argument('--out', action='store', help="Output .fasta file name (or stdout if not set")
TfamOptionsParser.add_argument('--action', required = True, action='store', help="Mutator function applied to the data")
TfamOptionsParser.add_argument('--list', action='store', help="List of entry names to use in filtering functions")

'''
    MAIN SCRIPT
'''
if __name__ == '__main__':

    run_args = TfamOptionsParser.parse_args(sys.argv[1:])

    input_source = sys.stdin
    if run_args.fasta is not None:
        if os.path.isfile(run_args.fasta):
            input_source = open(run_args.fasta, 'r')

    fasta_data = fasta_reader(input_source)
    filter_dict = {}
    if run_args.list is not None:
        if os.path.isfile(run_args.list):
            filter_dict = load_filter_dict(run_args.list)

    output_dest = sys.stdout
    output_is_file = False
    if run_args.out is not None:
        output_dest = open(run_args.out, 'w')
        output_is_file = True

    if run_args.action is None:
        printer('Error! Argument "--action" is mandatory\n', message_type = 'marked', stream = sys.stderr)
        exit(1)

    options = {
        'filter_dict' : filter_dict
    }

    action_func = ACTIONS[run_args.action]
    action_func(fasta_data, **options)
    write_fasta(fasta_data, output_dest)
    if output_is_file:
        output_dest.close()
