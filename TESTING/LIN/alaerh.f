      void alaerh(PATH, SUBNAM, INFO, INFOE, OPTS, M, N, KL, KU, N5, IMAT, NFAIL, NERRS, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      List<String>       SUBNAM;
      List<String>       OPTS;
      int                IMAT, INFO, INFOE, KL, KU, M, N, N5, NERRS, NFAIL, NOUT;
      // ..

// =====================================================================

      // .. Local Scalars ..
      String             UPLO;
      String             P2;
      String             C3;
      // ..
      // .. External Functions ..
      bool               LSAME, LSAMEN;
      // EXTERNAL LSAME, LSAMEN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN_TRIM
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAHD
      // ..
      // .. Executable Statements ..

      if (INFO == 0) return;
      P2 = PATH( 2: 3 );
      C3 = SUBNAM( 4: 6 );

      // Print the header if this is the first error message.

      if ( NFAIL == 0 && NERRS == 0 ) {
         if ( LSAMEN( 3, C3, 'SV ' ) || LSAMEN( 3, C3, 'SVX' ) ) {
            aladhd(NOUT, PATH );
         } else {
            alahd(NOUT, PATH );
         }
      }
      NERRS = NERRS + 1;

      // Print the message detailing the error and form of recovery,
      // if any.

      if ( LSAMEN( 2, P2, 'GE' ) ) {

         // xGE:  General matrices

         if ( LSAMEN( 3, C3, 'TRF' ) ) {
            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9988 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, M, N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9975 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, N5, IMAT;
            }
            if (INFO != 0) WRITE( NOUT, FMT = 9949 );

         } else if ( LSAMEN( 3, C3, 'SV ' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9984 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9970 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9992 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9997 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'TRI' ) ) {

            WRITE( NOUT, FMT = 9971 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, N, N5, IMAT;

         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) ) {

            WRITE( NOUT, FMT = 9978 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, IMAT;

         } else if ( LSAMEN( 3, C3, 'CON' ) ) {

            WRITE( NOUT, FMT = 9969 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), M, IMAT;

         } else if ( LSAMEN( 3, C3, 'LS ' ) ) {

            WRITE( NOUT, FMT = 9965 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), M, N, KL, N5, IMAT;

         } else if ( LSAMEN( 3, C3, 'LSX' ) || LSAMEN( 3, C3, 'LSS' ) ) {

            WRITE( NOUT, FMT = 9974 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, KL, N5, IMAT;

         } else {

            WRITE( NOUT, FMT = 9963 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), M, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'GB' ) ) {

         // xGB:  General band matrices

         if ( LSAMEN( 3, C3, 'TRF' ) ) {
            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9989 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, M, N, KL, KU, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9976 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, KL, KU, N5, IMAT;
            }
            if (INFO != 0) WRITE( NOUT, FMT = 9949 );

         } else if ( LSAMEN( 3, C3, 'SV ' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9986 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, N, KL, KU, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9972 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, N, KL, KU, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9993 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, OPTS( 1: 1 ), OPTS( 2: 2 ), N, KL, KU, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9998 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), N, KL, KU, N5, IMAT;
            }

         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) ) {

            WRITE( NOUT, FMT = 9977 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, KL, KU, IMAT;

         } else if ( LSAMEN( 3, C3, 'CON' ) ) {

            WRITE( NOUT, FMT = 9968 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), M, KL, KU, IMAT;

         } else {

            WRITE( NOUT, FMT = 9964 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), M, KL, KU, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'GT' ) ) {

         // xGT:  General tridiagonal matrices

         if ( LSAMEN( 3, C3, 'TRF' ) ) {
            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9987 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, N, IMAT;
            } else {
               WRITE( NOUT, FMT = 9973 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, N, IMAT;
            }
            if (INFO != 0) WRITE( NOUT, FMT = 9949 );

         } else if ( LSAMEN( 3, C3, 'SV ' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9984 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9970 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9992 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9997 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'CON' ) ) {

            WRITE( NOUT, FMT = 9969 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), M, IMAT;

         } else {

            WRITE( NOUT, FMT = 9963 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), M, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'PO' ) ) {

         // xPO:  Symmetric or Hermitian positive definite matrices

         UPLO = OPTS( 1: 1 );
         if ( LSAMEN( 3, C3, 'TRF' ) ) {
            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9980 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, UPLO, M, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9956 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, N5, IMAT;
            }
            if (INFO != 0) WRITE( NOUT, FMT = 9949 );

         } else if ( LSAMEN( 3, C3, 'SV ' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9979 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, UPLO, N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9955 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9990 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9995 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'TRI' ) ) {

            WRITE( NOUT, FMT = 9956 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, N5, IMAT;

         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) || LSAMEN( 3, C3, 'CON' ) ) {

            WRITE( NOUT, FMT = 9960 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, IMAT;

         } else {

            WRITE( NOUT, FMT = 9955 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'PS' ) ) {

         // xPS:  Symmetric or Hermitian positive semi-definite matrices

         UPLO = OPTS( 1: 1 );
         if ( LSAMEN( 3, C3, 'TRF' ) ) {
            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9980 )SUBNAM, INFO, INFOE, UPLO, M, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9956 )SUBNAM, INFO, UPLO, M, N5, IMAT;
            }
            if (INFO != 0) WRITE( NOUT, FMT = 9949 );

         } else if ( LSAMEN( 3, C3, 'SV ' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9979 )SUBNAM, INFO, INFOE, UPLO, N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9955 )SUBNAM, INFO, UPLO, N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9990 )SUBNAM, INFO, INFOE, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9995 )SUBNAM, INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'TRI' ) ) {

            WRITE( NOUT, FMT = 9956 )SUBNAM, INFO, UPLO, M, N5, IMAT;

         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMT' ) || LSAMEN( 3, C3, 'CON' ) ) {

            WRITE( NOUT, FMT = 9960 )SUBNAM, INFO, UPLO, M, IMAT;

         } else {

            WRITE( NOUT, FMT = 9955 )SUBNAM, INFO, UPLO, M, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'SY' ) || LSAMEN( 2, P2, 'SR' ) || LSAMEN( 2, P2, 'SK' ) || LSAMEN( 2, P2, 'HE' ) || LSAMEN( 2, P2, 'HR' ) || LSAMEN( 2, P2, 'HK' ) || LSAMEN( 2, P2, 'HA' ) ) {

         // xSY: symmetric indefinite matrices
              // with partial (Bunch-Kaufman) pivoting;
         // xSR: symmetric indefinite matrices
              // with rook (bounded Bunch-Kaufman) pivoting;
         // xSK: symmetric indefinite matrices
              // with rook (bounded Bunch-Kaufman) pivoting,
              // new storage format;
         // xHE: Hermitian indefinite matrices
              // with partial (Bunch-Kaufman) pivoting.
         // xHR: Hermitian indefinite matrices
              // with rook (bounded Bunch-Kaufman) pivoting;
         // xHK: Hermitian indefinite matrices
              // with rook (bounded Bunch-Kaufman) pivoting,
              // new storage format;
         // xHA: Hermitian matrices
              // Aasen Algorithm

         UPLO = OPTS( 1: 1 );
         if ( LSAMEN( 3, C3, 'TRF' ) ) {
            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9980 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, UPLO, M, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9956 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, N5, IMAT;
            }
            if (INFO != 0) WRITE( NOUT, FMT = 9949 );

         } else if ( LSAMEN( 2, C3, 'SV' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9979 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, UPLO, N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9955 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9990 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9995 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            }

         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) || LSAMEN( 3, C3, 'TRI' ) || LSAMEN( 3, C3, 'CON' ) ) {

            WRITE( NOUT, FMT = 9960 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, IMAT;

         } else {

            WRITE( NOUT, FMT = 9955 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'PP' ) || LSAMEN( 2, P2, 'SP' ) || LSAMEN( 2, P2, 'HP' ) ) {

         // xPP, xHP, or xSP:  Symmetric or Hermitian packed matrices

         UPLO = OPTS( 1: 1 );
         if ( LSAMEN( 3, C3, 'TRF' ) ) {
            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9983 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, UPLO, M, IMAT;
            } else {
               WRITE( NOUT, FMT = 9960 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, IMAT;
            }
            if (INFO != 0) WRITE( NOUT, FMT = 9949 );

         } else if ( LSAMEN( 3, C3, 'SV ' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9979 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, UPLO, N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9955 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9990 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9995 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), N, N5, IMAT;
            }

         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) || LSAMEN( 3, C3, 'TRI' ) || LSAMEN( 3, C3, 'CON' ) ) {

            WRITE( NOUT, FMT = 9960 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, IMAT;

         } else {

            WRITE( NOUT, FMT = 9955 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'PB' ) ) {

         // xPB:  Symmetric (Hermitian) positive definite band matrix

         UPLO = OPTS( 1: 1 );
         if ( LSAMEN( 3, C3, 'TRF' ) ) {
            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9982 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, UPLO, M, KL, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9958 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, KL, N5, IMAT;
            }
            if (INFO != 0) WRITE( NOUT, FMT = 9949 );

         } else if ( LSAMEN( 3, C3, 'SV ' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9981 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, UPLO, N, KL, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9957 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, N, KL, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9991 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, OPTS( 1: 1 ), OPTS( 2: 2 ), N, KL, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9996 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), N, KL, N5, IMAT;
            }

         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) || LSAMEN( 3, C3, 'CON' ) ) {

            WRITE( NOUT, FMT = 9959 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, KL, IMAT;

         } else {

            WRITE( NOUT, FMT = 9957 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, UPLO, M, KL, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'PT' ) ) {

         // xPT:  Positive definite tridiagonal matrices

         if ( LSAMEN( 3, C3, 'TRF' ) ) {
            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9987 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, N, IMAT;
            } else {
               WRITE( NOUT, FMT = 9973 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, N, IMAT;
            }
            if (INFO != 0) WRITE( NOUT, FMT = 9949 );

         } else if ( LSAMEN( 3, C3, 'SV ' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9984 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9970 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'SVX' ) ) {

            if ( INFO != INFOE && INFOE != 0 ) {
               WRITE( NOUT, FMT = 9994 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, OPTS( 1: 1 ), N, N5, IMAT;
            } else {
               WRITE( NOUT, FMT = 9999 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), N, N5, IMAT;
            }

         } else if ( LSAMEN( 3, C3, 'CON' ) ) {

            if( LSAME( SUBNAM( 1: 1 ), 'S' ) || LSAME( SUBNAM( 1: 1 ), 'D' ) ) {
               WRITE( NOUT, FMT = 9973 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, IMAT;
            } else {
               WRITE( NOUT, FMT = 9969 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), M, IMAT;
            }

         } else {

            WRITE( NOUT, FMT = 9963 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), M, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'TR' ) ) {

         // xTR:  Triangular matrix

         if ( LSAMEN( 3, C3, 'TRI' ) ) {
            WRITE( NOUT, FMT = 9961 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), M, N5, IMAT;
         } else if ( LSAMEN( 3, C3, 'CON' ) ) {
            WRITE( NOUT, FMT = 9967 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), OPTS( 3: 3 ), M, IMAT;
         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATRS' ) ) {
            WRITE( NOUT, FMT = 9952 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), OPTS( 3: 3 ), OPTS( 4: 4 ), M, IMAT;
         } else {
            WRITE( NOUT, FMT = 9953 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), OPTS( 3: 3 ), M, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'TP' ) ) {

         // xTP:  Triangular packed matrix

         if ( LSAMEN( 3, C3, 'TRI' ) ) {
            WRITE( NOUT, FMT = 9962 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), M, IMAT;
         } else if ( LSAMEN( 3, C3, 'CON' ) ) {
            WRITE( NOUT, FMT = 9967 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), OPTS( 3: 3 ), M, IMAT;
         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATPS' ) ) {
            WRITE( NOUT, FMT = 9952 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), OPTS( 3: 3 ), OPTS( 4: 4 ), M, IMAT;
         } else {
            WRITE( NOUT, FMT = 9953 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), OPTS( 3: 3 ), M, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'TB' ) ) {

         // xTB:  Triangular band matrix

         if ( LSAMEN( 3, C3, 'CON' ) ) {
            WRITE( NOUT, FMT = 9966 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), OPTS( 3: 3 ), M, KL, IMAT;
         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATBS' ) ) {
            WRITE( NOUT, FMT = 9951 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), OPTS( 3: 3 ), OPTS( 4: 4 ), M, KL, IMAT;
         } else {
            WRITE( NOUT, FMT = 9954 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, OPTS( 1: 1 ), OPTS( 2: 2 ), OPTS( 3: 3 ), M, KL, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'QR' ) ) {

         // xQR:  QR factorization

         if ( LSAMEN( 3, C3, 'QRS' ) ) {
            WRITE( NOUT, FMT = 9974 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, KL, N5, IMAT;
         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) ) {
            WRITE( NOUT, FMT = 9978 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'QK' ) ) {

         // xQK:  truncated QR factorization with pivoting

         if ( LSAMEN( 7, SUBNAM( 2: 8 ), 'GEQP3RK' )  ) {
            WRITE( NOUT, FMT = 9930 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, KL, N5, IMAT;
         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) ) {
            WRITE( NOUT, FMT = 9978 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'LQ' ) ) {

         // xLQ:  LQ factorization

         if ( LSAMEN( 3, C3, 'LQS' ) ) {
            WRITE( NOUT, FMT = 9974 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, KL, N5, IMAT;
         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) ) {
            WRITE( NOUT, FMT = 9978 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'QL' ) ) {

         // xQL:  QL factorization

         if ( LSAMEN( 3, C3, 'QLS' ) ) {
            WRITE( NOUT, FMT = 9974 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, KL, N5, IMAT;
         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) ) {
            WRITE( NOUT, FMT = 9978 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'RQ' ) ) {

         // xRQ:  RQ factorization

         if ( LSAMEN( 3, C3, 'RQS' ) ) {
            WRITE( NOUT, FMT = 9974 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, KL, N5, IMAT;
         } else if ( LSAMEN( 5, SUBNAM( 2: 6 ), 'LATMS' ) ) {
            WRITE( NOUT, FMT = 9978 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'LU' ) ) {

         if ( INFO != INFOE && INFOE != 0 ) {
            WRITE( NOUT, FMT = 9988 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, M, N, N5, IMAT;
         } else {
            WRITE( NOUT, FMT = 9975 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N, N5, IMAT;
         }

      } else if ( LSAMEN( 2, P2, 'CH' ) ) {

         if ( INFO != INFOE && INFOE != 0 ) {
            WRITE( NOUT, FMT = 9985 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, INFOE, M, N5, IMAT;
         } else {
            WRITE( NOUT, FMT = 9971 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO, M, N5, IMAT;
         }

      } else {

         // Print a generic message if the path is unknown.

         WRITE( NOUT, FMT = 9950 ) SUBNAM(1:LEN_TRIM( SUBNAM )), INFO;
      }

      // Description of error message (alphabetical, left to right)

      // SUBNAM, INFO, FACT, N, NRHS, IMAT

 9999 FORMAT( ' *** Error code from ', A, '=', I5, ', FACT=''', A1, ''', N=', I5, ', NRHS=', I4, ', type ', I2 );

      // SUBNAM, INFO, FACT, TRANS, N, KL, KU, NRHS, IMAT

 9998 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5, ', KL=', I5, ', KU=', I5, ', NRHS=', I4, ', type ', I1 );

      // SUBNAM, INFO, FACT, TRANS, N, NRHS, IMAT

 9997 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> FACT=''', A1, ''', TRANS=''', A1, ''', N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, FACT, UPLO, N, KD, NRHS, IMAT

 9996 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, ', KD=', I5, ', NRHS=', I4, ', type ', I2 );

      // SUBNAM, INFO, FACT, UPLO, N, NRHS, IMAT

 9995 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> FACT=''', A1, ''', UPLO=''', A1, ''', N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, FACT, N, NRHS, IMAT

 9994 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> FACT=''', A1, ''', N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, FACT, TRANS, N, KL, KU, NRHS, IMAT

 9993 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> FACT=''', A1, ''', TRANS=''', A1, ''', N=', I5, ', KL=', I5, ', KU=', I5, ', NRHS=', I4, ', type ', I1 );

      // SUBNAM, INFO, INFOE, FACT, TRANS, N, NRHS, IMAT

 9992 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> FACT=''', A1, ''', TRANS=''', A1, ''', N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, FACT, UPLO, N, KD, NRHS, IMAT

 9991 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> FACT=''', A1, ''', UPLO=''', A1, ''', N=', I5, ', KD=', I5, ', NRHS=', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, FACT, UPLO, N, NRHS, IMAT

 9990 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> FACT=''', A1, ''', UPLO=''', A1, ''', N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, M, N, KL, KU, NB, IMAT

 9989 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> M = ', I5, ', N =', I5, ', KL =', I5, ', KU =', I5, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, M, N, NB, IMAT

 9988 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> M =', I5, ', N =', I5, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, N, IMAT

 9987 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, ' for N=', I5, ', type ', I2 );

      // SUBNAM, INFO, INFOE, N, KL, KU, NRHS, IMAT

 9986 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> N =', I5, ', KL =', I5, ', KU =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, N, NB, IMAT

 9985 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> N =', I5, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, N, NRHS, IMAT

 9984 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, UPLO, N, IMAT

 9983 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> UPLO = ''', A1, ''', N =', I5, ', type ', I2 );

      // SUBNAM, INFO, INFOE, UPLO, N, KD, NB, IMAT

 9982 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> UPLO = ''', A1, ''', N =', I5, ', KD =', I5, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, UPLO, N, KD, NRHS, IMAT

 9981 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> UPLO=''', A1, ''', N =', I5, ', KD =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, UPLO, N, NB, IMAT

 9980 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, INFOE, UPLO, N, NRHS, IMAT

 9979 FORMAT( ' *** ', A, ' returned with INFO =', I5, ' instead of ', I2, / ' ==> UPLO = ''', A1, ''', N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, M, N, IMAT

 9978 FORMAT( ' *** Error code from ', A, ' =', I5, ' for M =', I5, ', N =', I5, ', type ', I2 );

      // SUBNAM, INFO, M, N, KL, KU, IMAT

 9977 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> M = ', I5, ', N =', I5, ', KL =', I5, ', KU =', I5, ', type ', I2 );

      // SUBNAM, INFO, M, N, KL, KU, NB, IMAT

 9976 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> M = ', I5, ', N =', I5, ', KL =', I5, ', KU =', I5, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, M, N, NB, IMAT

 9975 FORMAT( ' *** Error code from ', A, '=', I5, ' for M=', I5, ', N=', I5, ', NB=', I4, ', type ', I2 );

      // SUBNAM, INFO, M, N, NRHS, NB, IMAT

 9974 FORMAT( ' *** Error code from ', A, '=', I5, / ' ==> M =', I5, ', N =', I5, ', NRHS =', I4, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, N, IMAT

 9973 FORMAT( ' *** Error code from ', A, ' =', I5, ' for N =', I5, ', type ', I2 );

      // SUBNAM, INFO, N, KL, KU, NRHS, IMAT

 9972 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> N =', I5, ', KL =', I5, ', KU =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, N, NB, IMAT

 9971 FORMAT( ' *** Error code from ', A, '=', I5, ' for N=', I5, ', NB=', I4, ', type ', I2 );

      // SUBNAM, INFO, N, NRHS, IMAT

 9970 FORMAT( ' *** Error code from ', A, ' =', I5, ' for N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, NORM, N, IMAT

 9969 FORMAT( ' *** Error code from ', A, ' =', I5, ' for NORM = ''', A1, ''', N =', I5, ', type ', I2 );

      // SUBNAM, INFO, NORM, N, KL, KU, IMAT

 9968 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> NORM =''', A1, ''', N =', I5, ', KL =', I5, ', KU =', I5, ', type ', I2 );

      // SUBNAM, INFO, NORM, UPLO, DIAG, N, IMAT

 9967 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> NORM=''', A1, ''', UPLO =''', A1, ''', DIAG=''', A1, ''', N =', I5, ', type ', I2 );

      // SUBNAM, INFO, NORM, UPLO, DIAG, N, KD, IMAT

 9966 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> NORM=''', A1, ''', UPLO =''', A1, ''', DIAG=''', A1, ''', N=', I5, ', KD=', I5, ', type ', I2 );

      // SUBNAM, INFO, TRANS, M, N, NRHS, NB, IMAT

 9965 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> TRANS = ''', A1, ''', M =', I5, ', N =', I5, ', NRHS =', I4, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, TRANS, N, KL, KU, NRHS, IMAT

 9964 FORMAT( ' *** Error code from ', A, '=', I5, / ' ==> TRANS=''', A1, ''', N =', I5, ', KL =', I5, ', KU =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, TRANS, N, NRHS, IMAT

 9963 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> TRANS = ''', A1, ''', N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, UPLO, DIAG, N, IMAT

 9962 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO=''', A1, ''', DIAG =''', A1, ''', N =', I5, ', type ', I2 );

      // SUBNAM, INFO, UPLO, DIAG, N, NB, IMAT

 9961 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO=''', A1, ''', DIAG =''', A1, ''', N =', I5, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, UPLO, N, IMAT

 9960 FORMAT( ' *** Error code from ', A, ' =', I5, ' for UPLO = ''', A1, ''', N =', I5, ', type ', I2 );

      // SUBNAM, INFO, UPLO, N, KD, IMAT

 9959 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO = ''', A1, ''', N =', I5, ', KD =', I5, ', type ', I2 );

      // SUBNAM, INFO, UPLO, N, KD, NB, IMAT

 9958 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO = ''', A1, ''', N =', I5, ', KD =', I5, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, UPLO, N, KD, NRHS, IMAT

 9957 FORMAT( ' *** Error code from ', A, '=', I5, / ' ==> UPLO = ''', A1, ''', N =', I5, ', KD =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, UPLO, N, NB, IMAT

 9956 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ', I2 );

      // SUBNAM, INFO, UPLO, N, NRHS, IMAT

 9955 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO = ''', A1, ''', N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, UPLO, TRANS, DIAG, N, KD, NRHS, IMAT

 9954 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1, ''', N=', I5, ', KD=', I5, ', NRHS=', I4, ', type ', I2 );

      // SUBNAM, INFO, UPLO, TRANS, DIAG, N, NRHS, IMAT

 9953 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1, ''', N =', I5, ', NRHS =', I4, ', type ', I2 );

      // SUBNAM, INFO, UPLO, TRANS, DIAG, NORMIN, N, IMAT

 9952 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1, ''', NORMIN=''', A1, ''', N =', I5, ', type ', I2 );

      // SUBNAM, INFO, UPLO, TRANS, DIAG, NORMIN, N, KD, IMAT

 9951 FORMAT( ' *** Error code from ', A, ' =', I5, / ' ==> UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1, ''', NORMIN=''', A1, ''', N=', I5, ', KD=', I5, ', type ', I2 );

      // Unknown type

 9950 FORMAT( ' *** Error code from ', A, ' =', I5 );

      // What we do next

 9949 FORMAT( ' ==> Doing only the condition estimate for this case' );

      // SUBNAM, INFO, M, N, NB, IMAT

 9930 FORMAT( ' *** Error code from ', A, '=', I5, / ' ==> M =', I5, ', N =', I5, ', NX =', I5, ', NB =', I4, ', type ', I2 );

      return;
      }
