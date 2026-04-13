import h5py
from deerlab import FitResult, UQResult, Model, __VERSION__
import os
import pathlib
import numpy as np
import json
import tomllib
import re

supported_formats = ['hdf5', 'json', 'toml']
supported_types = ['FitResult', 'UQResult', ]


#=======================================================================================
#                              Serialization helpers (JSON / TOML)
#=======================================================================================

def _to_serializable(obj):
    """Recursively convert Python/numpy objects to JSON/TOML-safe primitives."""
    if isinstance(obj, np.ndarray):
        if np.issubdtype(obj.dtype, np.complexfloating):
            return {'__ndarray__': True, 'real': obj.real.tolist(), 'imag': obj.imag.tolist(), 'dtype': str(obj.dtype)}
        return {'__ndarray__': True, 'data': obj.tolist(), 'dtype': str(obj.dtype)}
    elif isinstance(obj, complex):
        return {'__complex__': True, 'real': obj.real, 'imag': obj.imag}
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, dict):
        return {str(k): _to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_to_serializable(i) for i in obj]
    elif obj is None:
        return {'__none__': True}
    return obj


def _from_serializable(obj):
    """Recursively restore numpy objects from serialized form."""
    if isinstance(obj, dict):
        if obj.get('__ndarray__') is True:
            if 'imag' in obj:
                return np.array(obj['real'], dtype=float) + 1j * np.array(obj['imag'], dtype=float)
            return np.array(obj['data'], dtype=obj.get('dtype', 'float64'))
        elif obj.get('__complex__') is True:
            return complex(obj['real'], obj['imag'])
        elif obj.get('__none__') is True:
            return None
        return {k: _from_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_from_serializable(i) for i in obj]
    return obj


#=======================================================================================
#                                         HDF5
#=======================================================================================

def _create_h5_element(group: h5py.Group, key, value):
    if isinstance(value, (Model, FitResult, UQResult)):
        return
    elif value is None:
        group.create_dataset(key, data='None')
    elif isinstance(value, list) and len(value) > 0 and isinstance(value[0], list):
        subgroup = group.create_group(key)
        subgroup.attrs['__type__'] = 'list'
        for i, sublist in enumerate(value):
            if isinstance(sublist, list):
                subgroup2 = subgroup.create_group(f'{key}_{i}')
                subgroup2.attrs['__type__'] = 'list'
                for j, item in enumerate(sublist):
                    _create_h5_element(subgroup2, f'{key}_{i}_{j}', item)
            else:
                subgroup.create_dataset(f'{key}_{i}', data=sublist)
    elif isinstance(value, dict):
        subgroup = group.create_group(key)
        for subkey, subvalue in value.items():
            _create_h5_element(subgroup, subkey, subvalue)
    else:
        group.create_dataset(key, data=value)


def _read_hdf5(filename):
    def _decode_h5_value(x):
        if isinstance(x, (bytes, np.bytes_)):
            s = x.decode('utf-8')
            return None if s == 'None' else s
        if isinstance(x, np.ndarray):
            if x.dtype.kind == 'S':
                return np.char.decode(x, 'utf-8')
            if x.dtype.kind == 'O':
                return np.array(
                    [i.decode('utf-8') if isinstance(i, (bytes, np.bytes_)) else i for i in x],
                    dtype=object,
                )
        return x

    def _h5_to_python(item):
        if isinstance(item, h5py.Group):
            if item.attrs.get('__type__') == 'list':
                return [_h5_to_python(item[k]) for k in item.keys()]
            return {subkey: _h5_to_python(item[subkey]) for subkey in item.keys()}
        else:
            return _decode_h5_value(item[()])

    with h5py.File(filename, 'r') as file:
        if file.attrs.get('format') != 'deerlab':
            raise ValueError(f"Unsupported format '{file.attrs.get('format')}' in file. Only 'deerlab' format is supported.")
        object_class = file.attrs['object_class']
        if object_class not in supported_types:
            raise ValueError(f"Unsupported object type '{object_class}' in file. Supported object types are: {supported_types}")

        raw_dict = {key: _h5_to_python(file[key]) for key in file.keys()}
        raw_dict['object_class'] = object_class

    return raw_dict


#=======================================================================================
#                                         JSON
#=======================================================================================

class _NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            if np.issubdtype(obj.dtype, np.complexfloating):
                return {'__ndarray__': True, 'real': obj.real.tolist(), 'imag': obj.imag.tolist(), 'dtype': str(obj.dtype)}
            return {'__ndarray__': True, 'data': obj.tolist(), 'dtype': str(obj.dtype)}
        elif isinstance(obj, complex):
            return {'__complex__': True, 'real': obj.real, 'imag': obj.imag}
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)


def _json_object_hook(d):
    if d.get('__ndarray__') is True:
        if 'imag' in d:
            return np.array(d['real'], dtype=float) + 1j * np.array(d['imag'], dtype=float)
        return np.array(d['data'], dtype=d.get('dtype', 'float64'))
    elif d.get('__complex__') is True:
        return complex(d['real'], d['imag'])
    elif d.get('__none__') is True:
        return None
    return d


def _save_json(filename, data_dict, object_class):
    payload = {'__format__': 'deerlab', '__version__': __VERSION__, '__object_class__': object_class}
    payload.update(data_dict)
    if isinstance(filename, (str, os.PathLike)):
        with open(filename, 'w') as f:
            json.dump(payload, f, cls=_NumpyEncoder, indent=2)
    else:
        json.dump(payload, filename, cls=_NumpyEncoder, indent=2)


def _load_json(filename):
    if isinstance(filename, (str, os.PathLike)):
        with open(filename, 'r') as f:
            payload = json.load(f, object_hook=_json_object_hook)
    else:
        payload = json.load(filename, object_hook=_json_object_hook)
    if payload.get('__format__') != 'deerlab':
        raise ValueError("File does not appear to be a deerlab JSON file.")
    object_class = payload.pop('__object_class__')
    payload.pop('__format__', None)
    payload.pop('__version__', None)
    if object_class not in supported_types:
        raise ValueError(f"Unsupported object type '{object_class}' in file.")
    payload['object_class'] = object_class
    return payload


#=======================================================================================
#                                         TOML
#=======================================================================================

def _toml_key(key):
    """Return a valid TOML key, quoting if necessary."""
    if re.match(r'^[A-Za-z0-9_-]+$', key):
        return key
    escaped = key.replace('\\', '\\\\').replace('"', '\\"')
    return f'"{escaped}"'


def _toml_value(val):
    """Serialize a value to an inline TOML value string."""
    if isinstance(val, bool):
        return 'true' if val else 'false'
    elif isinstance(val, int):
        return str(val)
    elif isinstance(val, float):
        if val != val:
            return 'nan'
        elif val == float('inf'):
            return 'inf'
        elif val == float('-inf'):
            return '-inf'
        return repr(val)
    elif isinstance(val, str):
        result = []
        for ch in val:
            if ch == '\\':
                result.append('\\\\')
            elif ch == '"':
                result.append('\\"')
            elif ch == '\n':
                result.append('\\n')
            elif ch == '\r':
                result.append('\\r')
            elif ch == '\t':
                result.append('\\t')
            elif ord(ch) < 0x20 or ord(ch) == 0x7f:
                result.append(f'\\u{ord(ch):04x}')
            else:
                result.append(ch)
        return '"' + ''.join(result) + '"'
    elif isinstance(val, (list, tuple)):
        return '[' + ', '.join(_toml_value(i) for i in val) + ']'
    elif isinstance(val, dict):
        items = ', '.join(f'{_toml_key(k)} = {_toml_value(v)}' for k, v in val.items())
        return '{' + items + '}'
    else:
        raise TypeError(f"Cannot serialize type {type(val).__name__} to TOML")


def _dict_to_toml(d):
    """Serialize a dict to a TOML string. Top-level dicts become [sections]."""
    scalars = []
    sections = {}
    for key, val in d.items():
        if isinstance(val, dict):
            sections[key] = val
        else:
            scalars.append(f'{_toml_key(key)} = {_toml_value(val)}')

    lines = scalars[:]
    for key, val in sections.items():
        lines.append(f'\n[{_toml_key(key)}]')
        for subkey, subval in val.items():
            lines.append(f'{_toml_key(subkey)} = {_toml_value(subval)}')

    return '\n'.join(lines) + '\n'


def _save_toml(filename, data_dict, object_class):
    meta = {'__format__': 'deerlab', '__version__': __VERSION__, '__object_class__': object_class}
    serializable = _to_serializable(data_dict)
    payload = {**meta, **serializable}
    content = _dict_to_toml(payload)
    if isinstance(filename, (str, os.PathLike)):
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(content)
    else:
        filename.write(content.encode('utf-8'))


def _load_toml(filename):
    if isinstance(filename, (str, os.PathLike)):
        with open(filename, 'rb') as f:
            payload = tomllib.load(f)
    else:
        content = filename.read()
        if isinstance(content, str):
            content = content.encode('utf-8')
        payload = tomllib.loads(content.decode('utf-8'))
    if payload.get('__format__') != 'deerlab':
        raise ValueError("File does not appear to be a deerlab TOML file.")
    object_class = payload.pop('__object_class__')
    payload.pop('__format__', None)
    payload.pop('__version__', None)
    if object_class not in supported_types:
        raise ValueError(f"Unsupported object type '{object_class}' in file.")
    result = _from_serializable(payload)
    result['object_class'] = object_class
    return result


#=======================================================================================
#                                         Saving
#=======================================================================================

def save(filename, object, format=None):
    """
    Saves a deerlab object (FitResult, UQresult, Model) to file.

    Parameters
    ----------
    filename : str, bytes, os.PathLike, or file-like object
        The name of the file to which the data is saved, or a path-like object.
        If the file already exists, it will be overwritten. If the file does not exist, a new file will be created.
        A file-like buffer (e.g. BytesIO) can be used for HDF5 and TOML; a text buffer for JSON.
    object : dl.FitResult, dl.UQResult, dl.Model
         The deerlab object to be saved.
    format : str, optional
        The format in which to save the data. One of 'hdf5', 'json', 'toml'.
        If not provided, the format is inferred from the file extension.
    """
    if format is None and isinstance(filename, (str, os.PathLike)):
        ext = pathlib.Path(filename).suffix.lower()
        if ext in ['.hdf5', '.h5']:
            format = 'hdf5'
        elif ext == '.json':
            format = 'json'
        elif ext == '.toml':
            format = 'toml'

    if format not in supported_formats:
        raise ValueError(f"Unsupported format '{format}'. Supported formats are: {supported_formats}")

    object_class = type(object).__name__
    if object_class not in supported_types:
        raise ValueError(f"Unsupported object type '{object_class}'. Supported object types are: {supported_types}")

    if format in ('hdf5', 'h5'):
        with h5py.File(filename, 'w') as file:
            file.attrs['format'] = 'deerlab'
            file.attrs['version'] = __VERSION__
            file.attrs['object_class'] = object_class
            for key, value in object.to_dict().items():
                _create_h5_element(file, key, value)
    elif format == 'json':
        _save_json(filename, object.to_dict(), object_class)
    elif format == 'toml':
        _save_toml(filename, object.to_dict(), object_class)


def json_dumps(object):
    """Convert a deerlab object to a JSON string."""
    if not isinstance(object, (FitResult, UQResult, Model)):
        raise ValueError("Only FitResult, UQResult, and Model objects can be serialized to JSON.")
    return json.dumps({'__format__': 'deerlab', '__version__': __VERSION__, '__object_class__': type(object).__name__, **object.to_dict()}, cls=_NumpyEncoder, indent=2)

def json_loads(json_string):
    """Convert a JSON string back to a deerlab object."""
    payload = json.loads(json_string, object_hook=_json_object_hook)
    if payload.get('__format__') != 'deerlab':
        raise ValueError("String does not appear to be a deerlab JSON string.")
    object_class = payload.pop('__object_class__')
    payload.pop('__format__', None)
    payload.pop('__version__', None)
    if object_class == 'FitResult':
        return FitResult.from_dict(payload)
    elif object_class == 'UQResult':
        return UQResult.from_dict(payload)
    elif object_class == 'Model':
        return Model.from_dict(payload)
    else:
        raise ValueError(f"Unsupported object type '{object_class}' in JSON string.")
#=======================================================================================
#                                         Loading
#=======================================================================================

def load(filename, format=None):
    """
    Loads a deerlab object (FitResult, UQresult, Model) from file.

    Parameters
    ----------
    filename : str, bytes, os.PathLike, or file-like object
        The name of the file from which the data is loaded, or a path-like object. The file must exist.
    format : str, optional
        The format of the file. One of 'hdf5', 'json', 'toml'.
        If not provided, the format is inferred from the file extension.

    Returns
    -------
    object : dl.FitResult, dl.UQResult, dl.Model
        The deerlab object that was loaded from file.
    """
    if format is None and isinstance(filename, (str, os.PathLike)):
        ext = pathlib.Path(filename).suffix.lower()
        if ext in ['.hdf5', '.h5']:
            file_format = 'hdf5'
        elif ext == '.json':
            file_format = 'json'
        elif ext == '.toml':
            file_format = 'toml'
        else:
            raise ValueError("Could not identify the file format. Please specify the format explicitly.")
    elif format is not None:
        file_format = format.lower()
        if file_format not in supported_formats:
            raise ValueError(f"Unsupported format '{file_format}'. Supported formats are: {supported_formats}")
    else:
        raise ValueError("Could not identify the file format. Please specify the format explicitly.")

    if file_format in ('hdf5', 'h5'):
        dict_output = _read_hdf5(filename)
    elif file_format == 'json':
        dict_output = _load_json(filename)
    elif file_format == 'toml':
        dict_output = _load_toml(filename)

    if dict_output['object_class'] == 'FitResult':
        dict_output.pop('object_class')
        return FitResult.from_dict(dict_output)
    elif dict_output['object_class'] == 'UQResult':
        dict_output.pop('object_class')
        return UQResult.from_dict(dict_output)
    elif dict_output['object_class'] == 'Model':
        dict_output.pop('object_class')
        return Model.from_dict(dict_output)