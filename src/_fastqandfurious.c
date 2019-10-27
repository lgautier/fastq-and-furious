#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <string.h>

#define POS_HEAD_BEG 0
#define POS_HEAD_END 1
#define POS_SEQ_BEG 2
#define POS_SEQ_END 3
#define POS_QUAL_BEG 4
#define POS_QUAL_END 5

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
  char * headerbeg_adr = (char *)memchr(blob_char + cur_offset, '@', blob.len - cur_offset - 1);
  if (headerbeg_adr == NULL) {
    headerbeg_adr = blob_char -1 ;    
  }
  posarray[POS_HEAD_BEG] = (Py_ssize_t) (headerbeg_adr - blob_char);

  if ((posarray[POS_HEAD_BEG] > 0) && (blob_char[posarray[POS_HEAD_BEG] - 1] == '+')) {
    PyErr_SetString(PyExc_Exception, "Lost sync' in input data.");
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return NULL;
  }
  if (posarray[POS_HEAD_BEG] == -1) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(0L);
  }

  cur_offset = posarray[POS_HEAD_BEG]+1;
  char * headerend_adr = (char *)memchr(blob_char + cur_offset, '\n', blob.len - cur_offset - 1);
  if (headerend_adr == NULL) {
    headerend_adr = blob_char - 1;
  }
  posarray[POS_HEAD_END] = (Py_ssize_t) (headerend_adr - blob_char);
  if (posarray[POS_HEAD_END] == -1) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(1L);
  }

  /* sequence */
  posarray[POS_SEQ_BEG] = posarray[POS_HEAD_END] + 1;
  if (posarray[POS_SEQ_BEG] == (blob.len-1)) {
    posarray[POS_SEQ_BEG] = -1;
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(2L);
  }
  cur_offset = posarray[POS_SEQ_BEG]+1;
  char * seqend_adr = (char *)memchr(blob_char + cur_offset, '\n', blob.len - cur_offset - 1);
  if (seqend_adr == NULL) {
    seqend_adr = blob_char - 1;
  }
  posarray[POS_SEQ_END] = (Py_ssize_t) (seqend_adr - blob_char);
  if (posarray[POS_SEQ_END] == (blob.len-1)) {
    posarray[POS_SEQ_END] = -1;
  }
  if (posarray[POS_SEQ_END] == -1) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(2L);
  }

  if ((posarray[POS_SEQ_END] + 1) >= blob.len) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(4L);
  }

  if (blob_char[posarray[POS_SEQ_END] + 1] != '+') {
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

  if ((posarray[POS_SEQ_END] + 2) >= blob.len) {
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(4L);
  }

  /* quality */

  /* POS_SEQ_END corresponds to the end of the sequence, including the EOL character \n.
   * If a FASTQ file withou quality header, there will only be '+\n' after this. 
   */
  if (blob_char[posarray[POS_SEQ_END] + 2] == '\n') {
    posarray[POS_QUAL_BEG] = posarray[POS_SEQ_END]+3;
    //        qualbeg_i = seqend_i+3
  } else {
    /* The header can optionally be repeated in the separator for quality. */
    Py_ssize_t lenheader = posarray[POS_HEAD_END] - posarray[POS_HEAD_BEG] + 1;
    if ((posarray[POS_HEAD_END] + lenheader) >= blob.len) {
      /* The buffer ends in the middle of an entry. */
      PyBuffer_Release(&blob);
      PyBuffer_Release(&posbuffer);
      return PyLong_FromLong(4L);
    } else if (blob_char[posarray[POS_SEQ_END] + lenheader] != '\n') {
      PyErr_Format(
		   PyExc_ValueError,
		   "Invalid quality header - expected end of line"
		   "(sequence header at blob[%lu:%lu], sequence at blob[%lu:%lu], "
		   "quality header at blob[%lu:%lu]).",
		   posarray[POS_HEAD_BEG]+1, posarray[POS_HEAD_END],
		   posarray[POS_SEQ_BEG], posarray[POS_SEQ_END],
		   posarray[POS_SEQ_END]+1, posarray[POS_SEQ_END]+lenheader
		   );
      PyBuffer_Release(&blob);
      PyBuffer_Release(&posbuffer);
      return NULL;
    } else {
      posarray[POS_QUAL_BEG] = posarray[POS_SEQ_END] + lenheader + 1;
    } 
  }

  if (posarray[POS_QUAL_BEG] >= blob.len) { // or posarray[3] == -1 ?
    posarray[POS_QUAL_BEG] = -1;
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(4L);
  }

  posarray[POS_QUAL_END] = posarray[POS_QUAL_BEG] + (posarray[POS_SEQ_END] - posarray[POS_SEQ_BEG]);
  if (posarray[POS_QUAL_END] >= blob.len) {
    posarray[POS_QUAL_END] = -1;
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return PyLong_FromLong(5L);
  } else if (blob_char[posarray[POS_QUAL_END]] != '\n') {
    PyErr_Format(PyExc_ValueError, "The last character in the quality line is not a new line (blob[%lu:%lu]).",
		 posarray[POS_QUAL_BEG], posarray[POS_QUAL_END]);
    PyBuffer_Release(&blob);
    PyBuffer_Release(&posbuffer);
    return NULL;
  }

  PyBuffer_Release(&blob);
  PyBuffer_Release(&posbuffer);
  return PyLong_FromLong(6L);
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
