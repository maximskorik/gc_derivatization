import pytest

from gc_meox_tms.__main__ import main
from mock import ANY, patch
from os.path import exists, join
from random import sample


@pytest.fixture
def input_data_path(test_dir):
    """Return the directory of the currently running test script as string type."""
    return join(test_dir, 'data/aldehydes.txt')


@pytest.fixture(params=[['-f'], ['-t'], ['-f', '-t']])
def output_params(request, tmp_path):
    """Return a list of output parameters."""
    args = []
    flat_path, tsv_path = None, None
    for flag in request.param:
        if flag == '-f':
            args.append(flag)
            args.append(flat_path := join(tmp_path, 'flat.txt'))
        elif flag == '-t':
            args.append(flag)
            args.append(tsv_path := join(tmp_path, 'tsv.txt'))
    yield args, flat_path, tsv_path


def test_cli_finishes(input_data_path):
    """Test if the main function works."""
    args = [input_data_path]
    exit_code = main(args)

    assert exit_code == 0


def test_cli_writes_files(input_data_path, output_params, tmp_path):
    args = output_params[0]
    args.append(input_data_path)

    flat_path = output_params[1]
    tsv_path = output_params[2]
    exit_code = main(args)

    assert exit_code == 0
    assert flat_path is None or exists(flat_path)
    assert tsv_path is None or exists(tsv_path)


@pytest.mark.parametrize('keep', [True, False])
@patch('gc_meox_tms.__main__.write_flat')
def test_keep_flag(mock, input_data_path, tmp_path, keep):
    """Test if the main function works with -k flag."""
    flat_path = join(tmp_path, 'flat.txt')
    args = [input_data_path, '-f', flat_path]
    if keep:
        args.append('-k')

    main(args)

    mock.assert_called_with(flat_path, ANY, keep)


@pytest.mark.parametrize('num_workers', sample(range(1, 10), 3))
@patch('gc_meox_tms.__main__.ProcessPoolExecutor')
def test_ncpu_flag(mock, input_data_path, num_workers):
    """Test if the main function works with -n flag."""
    args = [input_data_path, '-n', str(num_workers)]
    main(args)

    mock.assert_called_with(max_workers=num_workers)


@pytest.mark.parametrize('repeats', sample(range(1, 50), 3))
@patch('gc_meox_tms.__main__.ProcessPoolExecutor.map')
def test_repeats_flag(mock, input_data_path, repeats):
    """Test if the main function works with -r flag."""
    args = [input_data_path, '-r', str(repeats)]
    main(args)

    for params in mock.call_args[0][1]:
        assert params[1] == repeats
