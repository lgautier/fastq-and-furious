#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <string.h>

PyDoc_STRVAR(entrypos_doc,
             "entrypos(blob, backlog, posbuffer) -> int\n\n"
             "Compute the positions for the next FASTQ entry given "
	     "- blob: a bytes-like object\n"
	     "- backlog: a bytes-like object\n"
	     "- posbuffer: a buffer able to store 6 positions");

static PyObject *
entrypos(PyObject * self, PyObject * args)
{
  Py_buffer blob;
  Py_ssize_t offset;
  Py_buffer posbuffer;

  if (!PyArg_ParseTuple(args, "s*Ly*", &blob, &offset, &posbuffer)) {
    return NULL;
  }

  const char * blob_char = (char *)blob.buf;
      
  if (posbuffer.itemsize != sizeof(signed long long)) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    PyErr_SetString(PyExc_ValueError, "The buffer must be of format type q.");
    return NULL;
  }

  Py_ssize_t cur_offset = offset;
  Py_ssize_t * posarray = (Py_ssize_t *) posbuffer.buf;

  /* initialize the buffer */
  for (Py_ssize_t i = 0; i < 6; i++) {
    posarray[i] = -1;
  }
  
  /* header */
  char * headerbeg_adr = (char *)memchr(blob_char + cur_offset, '@', blob.len - cur_offset - 1);
  if (headerbeg_adr == NULL) {
    headerbeg_adr = blob_char -1 ;    
  }
  posarray[0] = (Py_ssize_t) (headerbeg_adr - blob_char);

  if ((posarray[0] > 0) && (blob_char[posarray[0] - 1] == '+')) {
    PyErr_SetString(PyExc_Exception, "Lost sync' in input data.");
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return NULL;
  }
  if (posarray[0] == -1) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(0L);
  }

  cur_offset = posarray[0]+1;
  char * headerend_adr = (char *)memchr(blob_char + cur_offset, '\n', blob.len - cur_offset - 1);
  if (headerend_adr == NULL) {
    headerend_adr = blob_char - 1;
  }
  posarray[1] = (Py_ssize_t) (headerend_adr - blob_char);
  if (posarray[1] == -1) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(1L);
  }

  /* sequence */
  posarray[2] = posarray[1] + 1;
  if (posarray[2] == (blob.len-1)) {
    posarray[2] = -1;
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(2L);
  }
  cur_offset = posarray[2]+1;
  char * seqend_adr = (char *)memchr(blob_char + cur_offset, '\n', blob.len - cur_offset - 1);
  if (seqend_adr == NULL) {
    seqend_adr = blob_char - 1;
  }
  posarray[3] = (Py_ssize_t) (seqend_adr - blob_char);
  if (posarray[3] == (blob.len-1)) {
    posarray[3] = -1;
  }
  if (posarray[3] == -1) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(2L);
  }

  if (blob_char[posarray[3] + 1] != '+') {
    /* multi-line FASTQ :/ */

    /* //char errstr[80]; */
    /* //PyOS_snprintf(errstr, (size_t)80, "", ...) */
    /* /\* FIXME: clean this up *\/ */
    /* printf("\noffset: %i\n", offset); */
    /* printf("posarray[2]: %i\n", posarray[2]); */
    /* printf("posarray[3]: %i\n", posarray[3]); */
    /* printf("header:\n"); */
    /* for(int i=posarray[0]; i<posarray[1]; i++)  */
    /*   printf("%c",blob_char[i]); */
    /* printf("\n---\n"); */
    /* printf("sequence:\n"); */
    /* for(int i=posarray[2]; i<posarray[3]; i++)  */
    /*   printf("%c",blob_char[i]); */
    /* printf("\n---\n"); */
    /* /\* --- *\/ */
    
    PyErr_SetString(PyExc_ValueError, "Multi-line FASTQ. Bye.");
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return NULL;
  }
  if ((posarray[3] + 2) >= blob.len) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(4L);
  }
  if (blob_char[posarray[3] + 2] != '\n') {
    PyErr_SetString(PyExc_ValueError,
		    "The character after the delimiter for quality should be a new line (the variant where the ID is repeated is not handled).");
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return NULL;
  }

  /* quality */
  posarray[4] = posarray[3]+3;
  if (posarray[4] == blob.len) { // or posarray[3] == -1 ?
    posarray[4] = -1;
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(4L);
  }

  posarray[5] = posarray[4] + (posarray[3] - posarray[2]);
  if (posarray[5] >= blob.len) {
    posarray[5] = -1;
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(5L);
  } else if (blob_char[posarray[5]] != '\n') {
    PyErr_SetString(PyExc_ValueError, "The last characted in the quality line is not a new line.");
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return NULL;
  }

  PyBuffer_Release(&blob);
  PyBuffer_Release(&posbuffer);
  return PyLong_FromLong(6L);
}

PyDoc_STRVAR(byteoffset_doc,
             "byteoffset(buf, offset)"
             "Offset the integer corresponding to each byte in the buffer **in-place**." 
	     "- buf: a bytes-like object\n"
	     "- offset: an integer\n");

static PyObject *
byteoffset(PyObject * self, PyObject * args)
{
  Py_buffer buf;
  signed char offset;

  if (!PyArg_ParseTuple(args, "y*h", &buf, &offset)) {
    return NULL;
  }

  if (buf.itemsize != sizeof(signed char)) {
    PyBuffer_Release(&buf);
    PyErr_SetString(PyExc_ValueError, "The buffer must be of format type b.");
    return NULL;
  }
  
  signed char * bufarray = (signed char *) buf.buf;
  const Py_ssize_t buflen = buf.len / buf.itemsize;
  /* initialize the buffer */
  for (Py_ssize_t i = 0; i < buflen; i++) {
    bufarray[i] += offset;
  }

  Py_RETURN_NONE;
}


static PyMethodDef fastqandfuriousModuleMethods[] = {
    {
      "entrypos", (PyCFunction)entrypos,
        METH_VARARGS, entrypos_doc,
    },
    {
      "byteoffset", (PyCFunction)byteoffset,
        METH_VARARGS, byteoffset_doc,
    },
    { NULL} // sentinel
};

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_fastqandfurious",
  "Utilities to handle FASTQ data.",
  -1,
  fastqandfuriousModuleMethods};

PyMODINIT_FUNC
PyInit__fastqandfurious(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);
    
    if (m == NULL) {
        return NULL;
    }

    return m;
}
