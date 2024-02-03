      SUBROUTINE ZDRVSG2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, LDB, D, D2, Z, LDZ, AB, BB, AP, BP, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, RESULT, INFO )

      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LDZ, LIWORK, LRWORK, NOUNIT, NSIZES, NTYPES, NWORK;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double             D( * ), D2( * ), RESULT( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AB( LDA, * ), AP( * ), B( LDB, * ), BB( LDB, * ), BP( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 10.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             UPLO;
      int                I, IBTYPE, IBUPLO, IINFO, IJ, IL, IMODE, ITEMP, ITYPE, IU, J, JCOL, JSIZE, JTYPE, KA, KA9, KB, KB9, M, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU, TEMP1, TEMP2;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLARND;
      // EXTERNAL LSAME, DLAMCH, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAFTS, DLASUM, XERBLA, ZHBGV, ZHBGVD, ZHBGVX, ZHEGV, ZHEGVD, ZHEGVX, ZHPGV, ZHPGVD, ZHPGVX, ZLACPY, ZLASET, ZLATMR, ZLATMS, ZSGT01, ZHEGV_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8, 6*9 /
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 6*1 /       DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 6*4 /
      // ..
      // .. Executable Statements ..

      // 1)      Check for errors

      NTESTT = 0
      INFO = 0

      BADNN = .FALSE.
      NMAX = 0
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE

      // Check for errors

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NTYPES.LT.0 ) {
         INFO = -3
      } else if ( LDA.LE.1 .OR. LDA.LT.NMAX ) {
         INFO = -9
      } else if ( LDZ.LE.1 .OR. LDZ.LT.NMAX ) {
         INFO = -16
      } else if ( 2*MAX( NMAX, 2 )**2.GT.NWORK ) {
         INFO = -21
      } else if ( 2*MAX( NMAX, 2 )**2.GT.LRWORK ) {
         INFO = -23
      } else if ( 2*MAX( NMAX, 2 )**2.GT.LIWORK ) {
         INFO = -25
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZDRVSG2STG', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

      // More Important constants

      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )

      DO 20 I = 1, 4
         ISEED2( I ) = ISEED( I )
   20 CONTINUE

      // Loop over sizes, types

      NERRS = 0
      NMATS = 0

      DO 650 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         ANINV = ONE / DBLE( MAX( 1, N ) )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         KA9 = 0
         KB9 = 0
         DO 640 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 640
            NMATS = NMATS + 1
            NTEST = 0

            DO 30 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   30       CONTINUE

            // 2)      Compute "A"

                    // Control parameters:

                // KMAGN  KMODE        KTYPE
            // =1  O(1)   clustered 1  zero
            // =2  large  clustered 2  identity
            // =3  small  exponential  (none)
            // =4         arithmetic   diagonal, w/ eigenvalues
            // =5         random log   hermitian, w/ eigenvalues
            // =6         random       (none)
            // =7                      random diagonal
            // =8                      random hermitian
            // =9                      banded, w/ eigenvalues

            IF( MTYPES.GT.MAXTYP ) GO TO 90

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE )

   40       CONTINUE
            ANORM = ONE
            GO TO 70

   50       CONTINUE
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 70

   60       CONTINUE
            ANORM = RTUNFL*N*ULPINV
            GO TO 70

   70       CONTINUE

            IINFO = 0
            COND = ULPINV

            // Special Matrices -- Identity & Jordan block

            if ( ITYPE.EQ.1 ) {

               // Zero

               KA = 0
               KB = 0
               CALL ZLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )

            } else if ( ITYPE.EQ.2 ) {

               // Identity

               KA = 0
               KB = 0
               CALL ZLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )
               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE

            } else if ( ITYPE.EQ.4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               KA = 0
               KB = 0
               CALL ZLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO )

            } else if ( ITYPE.EQ.5 ) {

               // Hermitian, eigenvalues specified

               KA = MAX( 0, N-1 )
               KB = KA
               CALL ZLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO )

            } else if ( ITYPE.EQ.7 ) {

               // Diagonal, random eigenvalues

               KA = 0
               KB = 0
               CALL ZLATMR( N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.8 ) {

               // Hermitian, random eigenvalues

               KA = MAX( 0, N-1 )
               KB = KA
               CALL ZLATMR( N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.9 ) {

               // Hermitian banded, eigenvalues specified

               // The following values are used for the half-bandwidths:

                 // ka = 1   kb = 1
                 // ka = 2   kb = 1
                 // ka = 2   kb = 2
                 // ka = 3   kb = 1
                 // ka = 3   kb = 2
                 // ka = 3   kb = 3

               KB9 = KB9 + 1
               if ( KB9.GT.KA9 ) {
                  KA9 = KA9 + 1
                  KB9 = 1
               }
               KA = MAX( 0, MIN( N-1, KA9 ) )
               KB = MAX( 0, MIN( N-1, KB9 ) )
               CALL ZLATMS( N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, KA, KA, 'N', A, LDA, WORK, IINFO )

            } else {

               IINFO = 1
            }

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

   90       CONTINUE

            ABSTOL = UNFL + UNFL
            if ( N.LE.1 ) {
               IL = 1
               IU = N
            } else {
               IL = 1 + INT( ( N-1 )*DLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( N-1 )*DLARND( 1, ISEED2 ) )
               if ( IL.GT.IU ) {
                  ITEMP = IL
                  IL = IU
                  IU = ITEMP
               }
            }

            // 3) Call ZHEGV, ZHPGV, ZHBGV, CHEGVD, CHPGVD, CHBGVD,
               // ZHEGVX, ZHPGVX and ZHBGVX, do tests.

            // loop over the three generalized problems
                  // IBTYPE = 1: A*x = (lambda)*B*x
                  // IBTYPE = 2: A*B*x = (lambda)*x
                  // IBTYPE = 3: B*A*x = (lambda)*x

            DO 630 IBTYPE = 1, 3

               // loop over the setting UPLO

               DO 620 IBUPLO = 1, 2
                  IF( IBUPLO.EQ.1 ) UPLO = 'U'                   IF( IBUPLO.EQ.2 ) UPLO = 'L'

                  // Generate random well-conditioned positive definite
                  // matrix B, of bandwidth not greater than that of A.

                  CALL ZLATMS( N, N, 'U', ISEED, 'P', RWORK, 5, TEN, ONE, KB, KB, UPLO, B, LDB, WORK( N+1 ), IINFO )

                  // Test ZHEGV

                  NTEST = NTEST + 1

                  CALL ZLACPY( ' ', N, N, A, LDA, Z, LDZ )
                  CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  CALL ZHEGV( IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, RWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  // Test ZHEGV_2STAGE

                  NTEST = NTEST + 1

                  CALL ZLACPY( ' ', N, N, A, LDA, Z, LDZ )
                  CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  CALL ZHEGV_2STAGE( IBTYPE, 'N', UPLO, N, Z, LDZ, BB, LDB, D2, WORK, NWORK, RWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 ) 'ZHEGV_2STAGE(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                   // CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z,
      // $                         LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  // Do Tests | D1 - D2 | / ( |D1| ulp )
                  // D1 computed using the standard 1-stage reduction as reference
                  // D2 computed using the 2-stage reduction

                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 151 J = 1, N
                     TEMP1 = MAX( TEMP1, ABS( D( J ) ),  ABS( D2( J ) ) )
                     TEMP2 = MAX( TEMP2, ABS( D( J )-D2( J ) ) )
  151             CONTINUE

                  RESULT( NTEST ) = TEMP2 /  MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

                  // Test ZHEGVD

                  NTEST = NTEST + 1

                  CALL ZLACPY( ' ', N, N, A, LDA, Z, LDZ )
                  CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  CALL ZHEGVD( IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  // Test ZHEGVX

                  NTEST = NTEST + 1

                  CALL ZLACPY( ' ', N, N, A, LDA, AB, LDA )
                  CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  CALL ZHEGVX( IBTYPE, 'V', 'A', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), IWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  NTEST = NTEST + 1

                  CALL ZLACPY( ' ', N, N, A, LDA, AB, LDA )
                  CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  // since we do not know the exact eigenvalues of this
                  // eigenpair, we just set VL and VU as constants.
                  // It is quite possible that there are no eigenvalues
                  // in this interval.

                  VL = ZERO
                  VU = ANORM
                  CALL ZHEGVX( IBTYPE, 'V', 'V', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), IWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  NTEST = NTEST + 1

                  CALL ZLACPY( ' ', N, N, A, LDA, AB, LDA )
                  CALL ZLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  CALL ZHEGVX( IBTYPE, 'V', 'I', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), IWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

  100             CONTINUE

                  // Test ZHPGV

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 120 J = 1, N
                        DO 110 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  110                   CONTINUE
  120                CONTINUE
                  } else {
                     IJ = 1
                     DO 140 J = 1, N
                        DO 130 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  130                   CONTINUE
  140                CONTINUE
                  }

                  CALL ZHPGV( IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, RWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  // Test ZHPGVD

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 160 J = 1, N
                        DO 150 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  150                   CONTINUE
  160                CONTINUE
                  } else {
                     IJ = 1
                     DO 180 J = 1, N
                        DO 170 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  170                   CONTINUE
  180                CONTINUE
                  }

                  CALL ZHPGVD( IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  // Test ZHPGVX

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 200 J = 1, N
                        DO 190 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  190                   CONTINUE
  200                CONTINUE
                  } else {
                     IJ = 1
                     DO 220 J = 1, N
                        DO 210 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  210                   CONTINUE
  220                CONTINUE
                  }

                  CALL ZHPGVX( IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, INFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 240 J = 1, N
                        DO 230 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  230                   CONTINUE
  240                CONTINUE
                  } else {
                     IJ = 1
                     DO 260 J = 1, N
                        DO 250 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  250                   CONTINUE
  260                CONTINUE
                  }

                  VL = ZERO
                  VU = ANORM
                  CALL ZHPGVX( IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, INFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 280 J = 1, N
                        DO 270 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  270                   CONTINUE
  280                CONTINUE
                  } else {
                     IJ = 1
                     DO 300 J = 1, N
                        DO 290 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  290                   CONTINUE
  300                CONTINUE
                  }

                  CALL ZHPGVX( IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, INFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

  310             CONTINUE

                  if ( IBTYPE.EQ.1 ) {

                     // TEST ZHBGV

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 340 J = 1, N
                           DO 320 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  320                      CONTINUE
                           DO 330 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  330                      CONTINUE
  340                   CONTINUE
                     } else {
                        DO 370 J = 1, N
                           DO 350 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  350                      CONTINUE
                           DO 360 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  360                      CONTINUE
  370                   CONTINUE
                     }

                     CALL ZHBGV( 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, RWORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                     // TEST ZHBGVD

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 400 J = 1, N
                           DO 380 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  380                      CONTINUE
                           DO 390 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  390                      CONTINUE
  400                   CONTINUE
                     } else {
                        DO 430 J = 1, N
                           DO 410 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  410                      CONTINUE
                           DO 420 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  420                      CONTINUE
  430                   CONTINUE
                     }

                     CALL ZHBGVD( 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                     // Test ZHBGVX

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 460 J = 1, N
                           DO 440 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  440                      CONTINUE
                           DO 450 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  450                      CONTINUE
  460                   CONTINUE
                     } else {
                        DO 490 J = 1, N
                           DO 470 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  470                      CONTINUE
                           DO 480 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  480                      CONTINUE
  490                   CONTINUE
                     }

                     CALL ZHBGVX( 'V', 'A', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 520 J = 1, N
                           DO 500 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  500                      CONTINUE
                           DO 510 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  510                      CONTINUE
  520                   CONTINUE
                     } else {
                        DO 550 J = 1, N
                           DO 530 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  530                      CONTINUE
                           DO 540 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  540                      CONTINUE
  550                   CONTINUE
                     }

                     VL = ZERO
                     VU = ANORM
                     CALL ZHBGVX( 'V', 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 580 J = 1, N
                           DO 560 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  560                      CONTINUE
                           DO 570 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  570                      CONTINUE
  580                   CONTINUE
                     } else {
                        DO 610 J = 1, N
                           DO 590 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  590                      CONTINUE
                           DO 600 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  600                      CONTINUE
  610                   CONTINUE
                     }

                     CALL ZHBGVX( 'V', 'I', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL ZSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  }

  620          CONTINUE
  630       CONTINUE

            // End of Loop -- Check for RESULT(j) > THRESH

            NTESTT = NTESTT + NTEST
            CALL DLAFTS( 'ZSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS )
  640    CONTINUE
  650 CONTINUE

      // Summary

      CALL DLASUM( 'ZSG', NOUNIT, NERRS, NTESTT )

      RETURN

 9999 FORMAT( ' ZDRVSG2STG: ', A, ' returned INFO=', I6, '.', / 9X,
     $  'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      // End of ZDRVSG2STG

      }
