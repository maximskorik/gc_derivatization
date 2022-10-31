import os
import pytest

from gc_meox_tms.__main__ import main
from mock import ANY, patch
from random import sample


@pytest.fixture(params=[['-f'], ['-t'], ['-f', '-t']])
def output_params(request, tmp_path):
    """Return a list of output parameters."""
    args = []
    flat_path, tsv_path = None, None
    for flag in request.param:
        if flag == '-f':
            args.append(flag)
            args.append(flat_path := os.path.join(tmp_path, 'flat.txt'))
        elif flag == '-t':
            args.append(flag)
            args.append(tsv_path := os.path.join(tmp_path, 'tsv.txt'))
    yield args, flat_path, tsv_path


def test_cli_finishes():
    """Test if the main function works."""
    args = ['data/aldehydes.txt']
    exit_code = main(args)

    assert exit_code == 0


def test_cli_writes_files(output_params, tmp_path):
    args = output_params[0]
    args.append('data/aldehydes.txt')

    flat_path = output_params[1]
    tsv_path = output_params[2]
    exit_code = main(args)

    assert exit_code == 0
    assert flat_path is None or os.path.exists(flat_path)
    assert tsv_path is None or os.path.exists(tsv_path)


@pytest.mark.parametrize('keep', [True, False])
@patch('gc_meox_tms.__main__.write_flat')
def test_keep_flag(mock, tmp_path, keep):
    """Test if the main function works with -k flag."""
    flat_path = os.path.join(tmp_path, 'flat.txt')
    args = ['data/alcohols.txt', '-f', flat_path]
    if keep:
        args.append('-k')

    main(args)

    mock.assert_called_with(flat_path, ANY, keep)


@pytest.mark.parametrize('num_workers', sample(range(1, 50), 3))
@patch('gc_meox_tms.__main__.ProcessPoolExecutor')
def test_ncpu_flag(mock, num_workers):
    """Test if the main function works with -n flag."""
    args = ['data/alcohols.txt', '-n', str(num_workers)]

    main(args)

    mock.assert_called_with(max_workers=num_workers)
