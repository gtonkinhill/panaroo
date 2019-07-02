import math


# set defaults based on mode if they have not been overwritten by the user
def set_default_args(args):

    n_samples = len(args.input_files)

    if args.mode == 'strict':
        if args.id is None:
            args.id = 0.95
        if args.family_threshold is None:
            args.family_threshold = 0.7
        if args.len_dif_percent is None:
            args.len_dif_percent = 0.95
        if args.min_trailing_support is None:
            args.min_trailing_support = max(2, math.ceil(0.05 * n_samples))
        if args.trailing_recursive is None:
            args.trailing_recursive = 99999999
        if args.min_edge_support_sv is None:
            args.min_edge_support_sv = max(2, math.ceil(0.01 * n_samples))
        if args.edge_support_diff is None:
            args.edge_support_diff = 0.05
        if args.remove_by_consensus is None:
            args.remove_by_consensus = True

    elif args.mode == 'moderate':
        if args.id is None:
            args.id = 0.95
        if args.family_threshold is None:
            args.family_threshold = 0.7
        if args.len_dif_percent is None:
            args.len_dif_percent = 0.95
        if args.min_trailing_support is None:
            args.min_trailing_support = max(2, math.ceil(0.01 * n_samples))
        if args.trailing_recursive is None:
            args.trailing_recursive = 99999999
        if args.min_edge_support_sv is None:
            args.min_edge_support_sv = max(2, math.ceil(0.01 * n_samples))
        if args.edge_support_diff is None:
            args.edge_support_diff = 0.01
        if args.remove_by_consensus is None:
            args.remove_by_consensus = False

    else:
        if args.id is None:
            args.id = 0.95
        if args.family_threshold is None:
            args.family_threshold = 0.7
        if args.len_dif_percent is None:
            args.len_dif_percent = 0.95
        if args.min_trailing_support is None:
            args.min_trailing_support = 2
        if args.trailing_recursive is None:
            args.trailing_recursive = 3
        if args.min_edge_support_sv is None:
            args.min_edge_support_sv = 2
        if args.edge_support_diff is None:
            args.edge_support_diff = 0.0
        if args.remove_by_consensus is None:
            args.remove_by_consensus = False

    return args
