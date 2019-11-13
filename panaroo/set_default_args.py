import math


# set defaults based on mode if they have not been overwritten by the user
def set_default_args(args):

    n_samples = len(args.input_files)

    if args.mode == 'strict':
        if args.id is None:
            args.id = 0.98
        if args.family_threshold is None:
            args.family_threshold = 0.7
        if args.len_dif_percent is None:
            args.len_dif_percent = 0.98
        if args.min_trailing_support is None:
            args.min_trailing_support = max(2, math.ceil(0.05 * n_samples))
        if args.trailing_recursive is None:
            args.trailing_recursive = 99999999
        if args.min_edge_support_sv is None:
            args.min_edge_support_sv = max(2, math.ceil(0.01 * n_samples))
        if args.remove_by_consensus is None:
            args.remove_by_consensus = True
        if args.edge_support_threshold is None:
            args.edge_support_threshold = max(2, math.ceil(0.01 * n_samples))

    elif args.mode == 'moderate':
        if args.id is None:
            args.id = 0.98
        if args.family_threshold is None:
            args.family_threshold = 0.7
        if args.len_dif_percent is None:
            args.len_dif_percent = 0.98
        if args.min_trailing_support is None:
            args.min_trailing_support = max(2, math.ceil(0.01 * n_samples))
        if args.trailing_recursive is None:
            args.trailing_recursive = 99999999
        if args.min_edge_support_sv is None:
            args.min_edge_support_sv = max(2, math.ceil(0.01 * n_samples))
        if args.remove_by_consensus is None:
            args.remove_by_consensus = False
        if args.edge_support_threshold is None:
            args.edge_support_threshold = max(2, math.ceil(0.01 * n_samples))

    else:
        if args.id is None:
            args.id = 0.98
        if args.family_threshold is None:
            args.family_threshold = 0.7
        if args.len_dif_percent is None:
            args.len_dif_percent = 0.98
        if args.min_trailing_support is None:
            args.min_trailing_support = 2
        if args.trailing_recursive is None:
            args.trailing_recursive = 0
        if args.min_edge_support_sv is None:
            args.min_edge_support_sv = 2
        if args.remove_by_consensus is None:
            args.remove_by_consensus = False
        if args.edge_support_threshold is None:
            args.edge_support_threshold = 0.0

    return args
