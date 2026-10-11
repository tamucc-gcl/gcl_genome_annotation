"""Run Terrier with the hierarchicalsoftmax 1.4.4 single-child inference fix.

Upstream skips softmax slots for single-child branches, but greedy_predictions
still indexes that missing slot when applying a threshold. Patch only that
condition in memory; keep upstream inference, model weights and CLI unchanged.
Source: rbturnbull/hierarchicalsoftmax, v1.4.4, inference.py.
"""
import inspect
import sys
from importlib.metadata import version


def install_single_child_fix():
    from hierarchicalsoftmax import inference

    if version('hierarchicalsoftmax') != '1.4.4':
        raise RuntimeError('Terrier compatibility fix requires hierarchicalsoftmax 1.4.4')
    original = inference.greedy_predictions
    source = inspect.getsource(original)
    condition = 'if threshold and predictions[node.softmax_start_index+prediction_child_index] < threshold:'
    if source.count(condition) != 1:
        raise RuntimeError('Unexpected upstream greedy_predictions source; review compatibility fix')
    source = source.replace(condition, condition.replace('if threshold', 'if len(node.children) > 1 and threshold'))
    namespace = dict(original.__globals__)
    exec(compile(source, '<hierarchicalsoftmax-single-child-fix>', 'exec'), namespace)
    inference.greedy_predictions = namespace['greedy_predictions']
    verify_single_child_fix(original, inference.greedy_predictions)
    print('Applied hierarchicalsoftmax 1.4.4 single-child threshold fix', file=sys.stderr)


def verify_single_child_fix(original, corrected):
    """Check the failing tree and confidence cutoff before loading model weights."""
    import torch
    from hierarchicalsoftmax import SoftmaxNode

    root = SoftmaxNode('root')
    first = SoftmaxNode('first', parent=root)
    second = SoftmaxNode('second', parent=root)
    leaf = SoftmaxNode('leaf', parent=second)
    root.set_indexes()
    # The single-child second branch starts just beyond the two score slots.
    confident = torch.tensor([[0.1, 0.9]])
    try:
        original(confident, root, threshold=0.7)
    except IndexError:
        pass
    else:
        raise RuntimeError('Expected upstream single-child regression did not reproduce')
    cases = [
        (confident, 0.7, None, leaf),
        (torch.tensor([[0.9, 0.1]]), 0.7, None, first),
        (torch.tensor([[0.4, 0.6]]), 0.7, None, root),
        (confident, 0.7, 1, second),
        (confident, None, None, leaf),
    ]
    for scores, threshold, depth, expected in cases:
        result = corrected(scores, root, threshold=threshold, max_depth=depth)
        if len(result) != 1 or result[0] is not expected:
            raise RuntimeError('Terrier single-child compatibility regression failed')


if __name__ == '__main__':
    install_single_child_fix()
    from terrier.apps import Terrier
    Terrier.main()
