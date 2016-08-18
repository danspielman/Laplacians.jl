const python = "/Users/serbanstan/anaconda/bin/python"
const libpython = "/Users/serbanstan/anaconda/lib/libpython2.7.dylib"
const pyprogramname = bytestring("/Users/serbanstan/anaconda/bin/python")
const pyversion_build = v"2.7.11"
const PYTHONHOME = bytestring("/Users/serbanstan/anaconda:/Users/serbanstan/anaconda")

"True if we are using the Python distribution in the Conda package."
const conda = false

const PyUnicode_AsUTF8String = :PyUnicodeUCS2_AsUTF8String
const PyUnicode_DecodeUTF8 = :PyUnicodeUCS2_DecodeUTF8

const PyString_FromStringAndSize = :PyString_FromStringAndSize
const PyString_AsStringAndSize = :PyString_AsStringAndSize
const PyString_Size = :PyString_Size
const PyString_Type = :PyString_Type
const PyInt_Type = :PyInt_Type
const PyInt_FromSize_t = :PyInt_FromSize_t
const PyInt_FromSsize_t = :PyInt_FromSsize_t
const PyInt_AsSsize_t = :PyInt_AsSsize_t

const Py_hash_t = Int64

const pyunicode_literals = false
