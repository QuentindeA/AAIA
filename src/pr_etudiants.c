/* PageRank
   The genetic.dat dataset comes from:
   http://www.cs.toronto.edu/~tsap/experiments/datasets/
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* allocate one object of given type */
#define NEW(type) ((type*)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type */
#define NEW_A(num,type) ((type*)calloc((size_t)(num),(size_t)sizeof(type)))

#define NBMULT 1000

#define ALPHA 0.0

typedef unsigned int u_int;

/* vector */
typedef struct
{
  u_int dim;
  double *e;
} VEC;

/* row of a sparse matrix */
typedef struct
{
  u_int  nnz;  /* # of non-zero (nz) value on this row */
  u_int  *col; /* column identifier for each nz value */
  double *val; /* value for each nz value */
} SROW;

/* sparse matrix */
typedef struct
{
  u_int m, n;
  SROW *row;
} SMAT;

/* v_get -- gets a VEC of dimension 'dim'
   Precondition: size >= 0
   Postcondition: initialized to zero */
VEC *v_get( u_int size )
{
  VEC *v;

  if( (v = NEW(VEC)) == (VEC *)NULL )
  {
    fprintf( stderr, "v_get memory error" );
    exit( -1 );
  }

  v->dim = size;
  if( (v->e = NEW_A(size,double)) == (double *)NULL )
  {
    free( v );
    fprintf( stderr, "v_get memory error" );
    exit( -1 );
  }

  return v;
}

/* v_free -- returns VEC & associated memory back to memory heap */
int v_free( VEC *v )
{
  if( v == (VEC *)NULL )
    return -1;

  if( v->e == (double *)NULL )
  {
    free( v );
  }
  else
  {
    free( v->e );
    free( v );
  }

  return 0;
}

/* sm_get -- gets an mxn sparse matrix by dynamic memory allocation
   Precondition: m>=0 && n>=0
   Postcondition: each row is empty*/
SMAT *sm_get( u_int m, u_int n )
{
  SMAT *M;
  u_int i;

  if( (M = NEW(SMAT)) == (SMAT *)NULL )
  {
    fprintf( stderr, "sm_get memory error" );
    exit( -1 );
  }

  M->m = m ; M->n = n;

  if( (M->row = NEW_A(m,SROW)) == (SROW *)NULL )
  {
    free( M );
    fprintf( stderr, "sm_get memory error" );
    exit( -1 );
  }

  for( i = 0 ; i < m ; i++ )
  {
    (M->row[i]).nnz = 0;
    (M->row[i]).col = (u_int *) NULL;
    (M->row[i]).val = (double *) NULL;
  }

  return M;
}

/* sm_free -- returns SMAT & associated memory back to memory heap */
int sm_free( SMAT *M )
{
  u_int i;
  SROW *ri;

  if( M == (SMAT *)NULL ) return -1;

  if( M->row == (SROW *)NULL )
  {
    free( M );
  }
  else
  {
    for( i = 0 ; i < M->m ; i++ )
    {
      ri = &(M->row[i]);
      if( ri->nnz > 0 )
      {
        free( ri->col );
        free( ri->val );
      }
    }
    free( M->row );
    free( M );
  }

  return 0;
}

/* sm_input -- file input of sparse matrix
   Precondition: will only work with a binary matrix. */
SMAT *sm_input( FILE *fp )
{
  SMAT *M;
  u_int *col; /* temp array to store the nz col of the current row */
  u_int r;    /* index of current row */
  int c;      /* index of current column */
  SROW *ri;   /* pointer to the current row in M */
  u_int m,n,i,j,k;

  /* get dimension */
  if( fscanf( fp, " SparseMatrix: %u by %u", &m, &n ) < 2 )
  {
    fprintf( stderr, "sm_input error reading dimensions" );
    exit( -1 );
  }

  if( (col = NEW_A(n,u_int)) == (u_int *)NULL )
  {
    fprintf( stderr, "sm_input memory error" );
    exit( -1 );
  }

  M = sm_get( m, n );

  /* get entries */
  for( i=0 ; i<m ; i++ )
  {
    if( fscanf( fp, " row %u:", &r ) < 1 )
    {
      fprintf( stderr, "sm_input error reading line %u", i );
      exit( -1 );
    }
    ri = &(M->row[i]);
    j = 0;
    for( ; ; )
    {
      if( fscanf( fp, "%d", &c ) < 1 )
      {
        fprintf( stderr, "sm_input error reading line %u col x", i );
        exit( -1 );
      }
      if( c < 0 ) break;
      col[j] = c;
      j++;
    } /* j is the number of nz value in row i */

    if( ( (ri->col = NEW_A(j,u_int)) == (u_int *)NULL ) && ( j!=0 ) )
    {
      fprintf( stderr, "sm_input memory error" );
      exit( -1 );
    }
    if( ( (ri->val = NEW_A(j,double)) == (double *)NULL ) && ( j!=0 ) )
    {
      fprintf( stderr, "sm_input memory error" );
      exit( -1 );
    }

    ri->nnz = j;

    for( k = 0 ; k < j ; k++ )
    {
      ri->col[k] = col[k];
      ri->val[k] = 1.0;
    }
  }

  free( col );

  return M;
}

/* sm_output -- file output of sparse matrix
   Postcondition: the result is not a valid entry for sm_input,
     since it also works for a non-binary matrix. */
void sm_output( FILE *fp, SMAT *M )
{
  u_int i,j;
  SROW *ri;

  fprintf( fp, "SparseMatrix: %d by %d\n", M->m, M->n );

  for( i = 0 ; i < M->m ; i++ )
  {
    fprintf( fp, "row %u: ", i );
    ri = &( M->row[i] );
    for( j = 0 ; j < ri->nnz ; j++ )
    {
      fprintf( fp, "%u:%1.5g ", ri->col[j], ri->val[j] );
    }
    fprintf( fp, "-1\n" );
  }
}

/* v_output -- file output of vector */
void v_output( FILE *fp, VEC *v )
{
  u_int i;

  fprintf( fp, "Vector: %d\n", v->dim );
  for( i = 0 ; i < v->dim ; i++ ) fprintf( fp, "%1.5g ", v->e[i] );
  putc( '\n', fp );
}

void m_to_h(SMAT *M)
{
    u_int sizeM = M->m;
    for(int i = 0; i < sizeM; i++)
    {
        for(int j = 0; j< M->row[i].nnz; j++)
        {
            M->row[i].val[j] = 1.0/M->row[i].nnz;
        }
    }

    /*for(int i = 0; i < sizeM; i++)
    {

        if (M->row[i].nnz == 0)
        {

            if( ( (M->row[i].col = NEW_A(1,u_int)) == (u_int *)NULL ) )
            {
              fprintf( stderr, "sm_input memory error" );
              exit( -1 );
            }
            if( ( (M->row[i].val = NEW_A(1,double)) == (double *)NULL ) )
            {
              fprintf( stderr, "sm_input memory error" );
              exit( -1 );
            }

            M->row[i].col[0] = i;
            M->row[i].val[0] = 1.0;
            M->row[i].nnz = 1;
        }
    }*/
}



VEC *createVec(int m)
{
    VEC *myVec = v_get(m);
    for(int i=0;i<m;i++)
    {
        myVec->e[i] = 1.0/m;
    }
    return myVec;
}

VEC *multiply(VEC *vec, SMAT *H )
{
    VEC *newVec;
    newVec = v_get(vec->dim);

    double standard = 1.0/H->m;

    for(int j=0; j<H->m;j++)
    {
        /*if(H->row[j].nnz == 0)
        {

            for(int k=0;k<H->m;k++)
            {
                newVec->e[k] += standard*vec->e[j];
            }
        }*/
        for(int k=0;k<H->row[j].nnz;k++)
        {
            newVec->e[H->row[j].col[k]] += H->row[j].val[k]*vec->e[j];
        }
    }


    return newVec;
}

VEC *surfeur(VEC *vec, SMAT *H)
{
    VEC *rk_old = multiply(vec, H);
    VEC *rk_alphaise = v_get(H->m);
    VEC *a = v_get(H->m);
    double standard = 1.0/H->m;
    VEC *rk_new = v_get(H->m);
    double indicator = 0.0;

    for (size_t i = 0; i < H->m; i++)
    {
        rk_alphaise->e[i] = rk_old->e[i]*ALPHA;
        if(H->row[i].nnz == 0)
        {
            a->e[i] = 1.0;
        }
    }

    for (size_t i = 0; i < H->m; i++) {
        indicator += a->e[i]*vec->e[i]*ALPHA;
    }


    for (size_t i = 0; i < H->m; i++) {
        rk_new->e[i] = rk_alphaise->e[i] + (indicator + 1.0 - ALPHA)*standard;

    }

    v_free(rk_old);
    v_free(a);
    v_free(rk_alphaise);
    return rk_new;


}

int main()
{
    FILE *fp;
    SMAT *SM;

    fp = fopen( "exemple.dat", "r" );
    SM = sm_input( fp );
    fclose( fp );

    m_to_h(SM);
    VEC *vec = createVec(SM->m);

    v_output(stdout, vec);
    printf("\n");
    for(int i=0; i<NBMULT; i++)
    {
      vec = surfeur(vec, SM);
    }

    v_output(stdout, vec);
    printf("\n");

    sm_output( stdout, SM );

    v_free(vec);
    sm_free( SM );

    return 0;
}
