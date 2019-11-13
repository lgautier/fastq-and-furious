#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <string.h>

#define INVALID -1
#define POS_HEAD_BEG 0
#define POS_HEAD_END 1
#define POS_SEQ_BEG 2
#define POS_SEQ_END 3
#define POS_QUAL_BEG 4
#define POS_QUAL_END 5
#define COMPLETE 6
#define MISSING_QUALHEADER_END 7

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

  char * blob_char = (char *)blob.buf;
      
  if (posbuffer.itemsize != sizeof(signed long long)) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    PyErr_SetString(PyExc_ValueError, "The buffer must be of format type q.");
    return NULL;
  }

  Py_ssize_t cur_offset = offset;
  /* posarray is an array with 6 offsets:
  * - 0: header, begin
  * - 1: header, end
  * - 2: sequence, begin
  * - 3: sequence, end
  * - 4: quality, begin
  * - 5: quality, end
  */
  Py_ssize_t * posarray = (Py_ssize_t *) posbuffer.buf;

  /* initialize the buffer */
  for (Py_ssize_t i = 0; i < 6; i++) {
    posarray[i] = -1;
  }
  
  /* header */  
  char * headerbeg_adr = (char *)memmem((void *)(blob_char + cur_offset), blob.len - cur_offset, (void *)("\n@"), 2);
  if (headerbeg_adr == NULL) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(POS_HEAD_BEG);
  }
  posarray[POS_HEAD_BEG] = (Py_ssize_t) (headerbeg_adr - blob_char + 1);

  cur_offset = posarray[POS_HEAD_BEG]+1;
  char * headerend_adr = (char *)memchr(blob_char + cur_offset, '\n', blob.len - cur_offset - 1);
  if (headerend_adr == NULL) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(POS_HEAD_END);
  }
  posarray[POS_HEAD_END] = (Py_ssize_t) (headerend_adr - blob_char);

  /* sequence */
  if (posarray[POS_HEAD_END]+1 >= blob.len) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(POS_SEQ_BEG);
  }
  posarray[POS_SEQ_BEG] = posarray[POS_HEAD_END] + 1;

  cur_offset = posarray[POS_SEQ_BEG]+1;
  char * seqend_adr = (char *)memmem((void *)(blob_char + cur_offset), blob.len - cur_offset, (void *)("\n+"), 2);
  if (seqend_adr == NULL) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(POS_SEQ_END);
  }
  posarray[POS_SEQ_END] = (Py_ssize_t) (seqend_adr - blob_char);

  /* Quality */
  if (posarray[POS_SEQ_END]+2 >= blob.len) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(MISSING_QUALHEADER_END);
  }
  cur_offset = posarray[POS_SEQ_END]+2;
  char * qualheadend_adr = (char *)memchr(blob_char + cur_offset, '\n', blob.len - cur_offset - 1);
  if (qualheadend_adr == NULL) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(MISSING_QUALHEADER_END);
  }
  if (
      ((qualheadend_adr - blob_char - posarray[POS_SEQ_END] - 1) > 1)
      &&
      ((qualheadend_adr - blob_char - posarray[POS_SEQ_END]) != (posarray[POS_HEAD_END] - posarray[POS_HEAD_BEG] + 1))
      ) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(INVALID);
  }
  
  Py_ssize_t qualbeg_i = qualheadend_adr - blob_char +1;
  
  if (qualbeg_i >= blob.len) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(POS_QUAL_BEG);
  } else {
    posarray[POS_QUAL_BEG] = qualbeg_i;
  }

  Py_ssize_t qualend_i = (posarray[POS_QUAL_BEG] + posarray[POS_SEQ_END] - posarray[POS_HEAD_END] - 1);
  if ((qualend_i+2) >= blob.len) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(POS_QUAL_END);
  } else {
    posarray[POS_QUAL_END] = qualend_i;
  }

  if (!(
	((qualend_i-posarray[POS_QUAL_BEG]) != (posarray[POS_SEQ_END]+1 - posarray[POS_HEAD_END]))
      ||
      ((! (qualend_i+2 >= blob.len)) && blob_char[qualend_i] == '\n' && blob_char[qualend_i+1] == '@')
      ||
      (blob.len - qualend_i)
	)
      ) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(INVALID);
  }
  PyBuffer_Release(&blob);
  PyBuffer_Release(&posbuffer);
  return PyLong_FromLong(COMPLETE);
}

PyDoc_STRVAR(arrayadd_b_doc,
             "arrayadd_b(a, offset)"
             "Add the integer value to each value in the array **in-place**." 
 	     "- a: an array (or any object with the buffer interface) with elements of type 'b'\n"
	     "- value: a signed integer\n");

static PyObject *
arrayadd_b(PyObject * self, PyObject * args)
{
  Py_buffer buf;
  signed char value;

  if (!PyArg_ParseTuple(args, "y*h", &buf, &value)) {
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
    bufarray[i] += value;
  }

  Py_RETURN_NONE;
}

PyDoc_STRVAR(arrayadd_q_doc,
             "arrayadd_q(a, offset)"
	     "Add the integer value to each value in the array **in-place**." 
 	     "- a: an array (or any object with the buffer interface) with elements of type 'q'\n"
	     "- value: a signed integer\n");

static PyObject *
arrayadd_q(PyObject * self, PyObject * args)
{
  Py_buffer buf;
  long long value;

  if (!PyArg_ParseTuple(args, "y*L", &buf, &value)) {
    return NULL;
  }

  if (buf.itemsize != sizeof(unsigned long long)) {
    PyBuffer_Release(&buf);
    PyErr_SetString(PyExc_ValueError, "The buffer must be of format type q.");
    return NULL;
  }
  
  unsigned long long * bufarray = (unsigned long long *) buf.buf;
  const Py_ssize_t buflen = buf.len / buf.itemsize;
  /* initialize the buffer */
  for (Py_ssize_t i = 0; i < buflen; i++) {
    bufarray[i] += value;
  }

  Py_RETURN_NONE;
}


static PyMethodDef fastqandfuriousModuleMethods[] = {
    {
      "entrypos", (PyCFunction)entrypos,
        METH_VARARGS, entrypos_doc,
    },
    {
      "arrayadd_b", (PyCFunction)arrayadd_b,
        METH_VARARGS, arrayadd_b_doc,
    },
    {
      "arrayadd_q", (PyCFunction)arrayadd_q,
        METH_VARARGS, arrayadd_q_doc,
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
