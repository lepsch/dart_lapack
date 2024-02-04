      void schkbd(NSIZES, MVAL, NVAL, NTYPES, DOTYPE, NRHS, ISEED, THRESH, A, LDA, BD, BE, S1, S2, X, LDX, Y, Z, Q, LDQ, PT, LDPT, U, VT, WORK, LWORK, IWORK, NOUT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDPT, LDQ, LDX, LWORK, NOUT, NRHS, NSIZES, NTYPES;
      double               THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * );
      double               A( LDA, * ), BD( * ), BE( * ), PT( LDPT, * ), Q( LDQ, * ), S1( * ), S2( * ), U( LDPT, * ), VT( LDPT, * ), WORK( * ), X( LDX, * ), Y( LDX, * ), Z( LDX, * );
      // ..

// ======================================================================

      // .. Parameters ..
      double               ZERO, ONE, TWO, HALF;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 0.5 ;
      int                MAXTYP;
      const              MAXTYP = 16 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN, BIDIAG;
      String             UPLO;
      String             PATH;
      int                I, IINFO, IL, IMODE, ITEMP, ITYPE, IU, IWBD, IWBE, IWBS, IWBZ, IWWORK, J, JCOL, JSIZE, JTYPE, LOG2UI, M, MINWRK, MMAX, MNMAX, MNMIN, MNMIN2, MQ, MTYPES, N, NFAIL, NMAX, NS1, NS2, NTEST;
      double               ABSTOL, AMNINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL, VL, VU;
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 ), IOLDSD( 4 ), ISEED2( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      double               DUM( 1 ), DUMMA( 1 ), RESULT( 40 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLARND, SSXT1;
      // EXTERNAL SLAMCH, SLARND, SSXT1
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASUM, SBDSDC, SBDSQR, SBDSVDX, SBDT01, SBDT02, SBDT03, SBDT04, SCOPY, SGEBRD, SGEMM, SLACPY, SLAHD2, SLASET, SLATMR, SLATMS, SORGBR, SORT01, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, EXP, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 9, 9, 9, 10 ];
      const KMAGN = [ 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 0 ];
      const KMODE = [ 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 0 ];
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0;

      BADMM = false;
      BADNN = false;
      MMAX = 1;
      NMAX = 1;
      MNMAX = 1;
      MINWRK = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         MMAX = max( MMAX, MVAL( J ) );
         if( MVAL( J ) < 0 ) BADMM = true;
         NMAX = max( NMAX, NVAL( J ) );
         if( NVAL( J ) < 0 ) BADNN = true;
         MNMAX = max( MNMAX, min( MVAL( J ), NVAL( J ) ) );
         MINWRK = max( MINWRK, 3*( MVAL( J )+NVAL( J ) ), MVAL( J )*( MVAL( J )+max( MVAL( J ), NVAL( J ), NRHS )+1 )+NVAL( J )*min( NVAL( J ), MVAL( J ) ) );
      } // 10

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1;
      } else if ( BADMM ) {
         INFO = -2;
      } else if ( BADNN ) {
         INFO = -3;
      } else if ( NTYPES < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -6;
      } else if ( LDA < MMAX ) {
         INFO = -11;
      } else if ( LDX < MMAX ) {
         INFO = -17;
      } else if ( LDQ < MMAX ) {
         INFO = -21;
      } else if ( LDPT < MNMAX ) {
         INFO = -23;
      } else if ( MINWRK > LWORK ) {
         INFO = -27;
      }

      if ( INFO != 0 ) {
         xerbla('SCHKBD', -INFO );
         return;
      }

      // Initialize constants

      PATH[1: 1] = 'Single precision';
      PATH[2: 3] = 'BD';
      NFAIL = 0;
      NTEST = 0;
      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = SLAMCH( 'Overflow' );
      ULP = SLAMCH( 'Precision' );
      ULPINV = ONE / ULP;
      LOG2UI = INT( LOG( ULPINV ) / LOG( TWO ) );
      RTUNFL = sqrt( UNFL );
      RTOVFL = sqrt( OVFL );
      INFOT = 0;
      ABSTOL = 2*UNFL;

      // Loop over sizes, types

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 300
         M = MVAL( JSIZE );
         N = NVAL( JSIZE );
         MNMIN = min( M, N );
         AMNINV = ONE / max( M, N, 1 );

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 290
            if( !DOTYPE( JTYPE ) ) GO TO 290;

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD[J] = ISEED( J );
            } // 20

            for (J = 1; J <= 34; J++) { // 30
               RESULT[J] = -ONE;
            } // 30

            UPLO = ' ';

            // Compute "A"

            // Control parameters:

            // KMAGN  KMODE        KTYPE
        // =1  O(1)   clustered 1  zero
        // =2  large  clustered 2  identity
        // =3  small  exponential  (none)
        // =4         arithmetic   diagonal, (w/ eigenvalues)
        // =5         random       symmetric, w/ eigenvalues
        // =6                      nonsymmetric, w/ singular values
        // =7                      random diagonal
        // =8                      random symmetric
        // =9                      random nonsymmetric
        // =10                     random bidiagonal (log. distrib.)

            if (MTYPES > MAXTYP) GO TO 100;

            ITYPE = KTYPE( JTYPE );
            IMODE = KMODE( JTYPE );

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE );

            } // 40
            ANORM = ONE;
            GO TO 70;

            } // 50
            ANORM = ( RTOVFL*ULP )*AMNINV;
            GO TO 70;

            } // 60
            ANORM = RTUNFL*max( M, N )*ULPINV;
            GO TO 70;

            } // 70

            slaset('Full', LDA, N, ZERO, ZERO, A, LDA );
            IINFO = 0;
            COND = ULPINV;

            BIDIAG = false;
            if ( ITYPE == 1 ) {

               // Zero matrix

               IINFO = 0;

            } else if ( ITYPE == 2 ) {

               // Identity

               for (JCOL = 1; JCOL <= MNMIN; JCOL++) { // 80
                  A[JCOL, JCOL] = ANORM;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               slatms(MNMIN, MNMIN, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( MNMIN+1 ), IINFO );

            } else if ( ITYPE == 5 ) {

               // Symmetric, eigenvalues specified

               slatms(MNMIN, MNMIN, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, M, N, 'N', A, LDA, WORK( MNMIN+1 ), IINFO );

            } else if ( ITYPE == 6 ) {

               // Nonsymmetric, singular values specified

               slatms(M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, M, N, 'N', A, LDA, WORK( MNMIN+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random entries

               slatmr(MNMIN, MNMIN, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( 2*MNMIN+1 ), 1, ONE, 'N', IWORK, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Symmetric, random entries

               slatmr(MNMIN, MNMIN, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( M+MNMIN+1 ), 1, ONE, 'N', IWORK, M, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // Nonsymmetric, random entries

               slatmr(M, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( M+MNMIN+1 ), 1, ONE, 'N', IWORK, M, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 10 ) {

               // Bidiagonal, random entries

               TEMP1 = -TWO*LOG( ULP );
               for (J = 1; J <= MNMIN; J++) { // 90
                  BD[J] = EXP( TEMP1*SLARND( 2, ISEED ) );
                  if (J < MNMIN) BE( J ) = EXP( TEMP1*SLARND( 2, ISEED ) );
               } // 90

               IINFO = 0;
               BIDIAG = true;
               if ( M >= N ) {
                  UPLO = 'U';
               } else {
                  UPLO = 'L';
               }
            } else {
               IINFO = 1;
            }

            if ( IINFO == 0 ) {

               // Generate Right-Hand Side

               if ( BIDIAG ) {
                  slatmr(MNMIN, NRHS, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( 2*MNMIN+1 ), 1, ONE, 'N', IWORK, MNMIN, NRHS, ZERO, ONE, 'NO', Y, LDX, IWORK, IINFO );
               } else {
                  slatmr(M, NRHS, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( M+1 ), 1, ONE, WORK( 2*M+1 ), 1, ONE, 'N', IWORK, M, NRHS, ZERO, ONE, 'NO', X, LDX, IWORK, IINFO );
               }
            }

            // Error Exit

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'Generator', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               return;
            }

            } // 100

            // Call SGEBRD and SORGBR to compute B, Q, and P, do tests.

            if ( !BIDIAG ) {

               // Compute transformations to reduce A to bidiagonal form:
               // B := Q' * A * P.

               slacpy(' ', M, N, A, LDA, Q, LDQ );
               sgebrd(M, N, Q, LDQ, BD, BE, WORK, WORK( MNMIN+1 ), WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO );

               // Check error code from SGEBRD.

               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9998 )'SGEBRD', IINFO, M, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               slacpy(' ', M, N, Q, LDQ, PT, LDPT );
               if ( M >= N ) {
                  UPLO = 'U';
               } else {
                  UPLO = 'L';
               }

               // Generate Q

               MQ = M;
               if (NRHS <= 0) MQ = MNMIN;
               sorgbr('Q', M, MQ, N, Q, LDQ, WORK, WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO );

               // Check error code from SORGBR.

               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9998 )'SORGBR(Q)', IINFO, M, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               // Generate P'

               sorgbr('P', MNMIN, N, M, PT, LDPT, WORK( MNMIN+1 ), WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO );

               // Check error code from SORGBR.

               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9998 )'SORGBR(P)', IINFO, M, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               // Apply Q' to an M by NRHS matrix X:  Y := Q' * X.

               sgemm('Transpose', 'No transpose', M, NRHS, M, ONE, Q, LDQ, X, LDX, ZERO, Y, LDX );

               // Test 1:  Check the decomposition A := Q * B * PT
                    // 2:  Check the orthogonality of Q
                    // 3:  Check the orthogonality of PT

               sbdt01(M, N, 1, A, LDA, Q, LDQ, BD, BE, PT, LDPT, WORK, RESULT( 1 ) );
               sort01('Columns', M, MQ, Q, LDQ, WORK, LWORK, RESULT( 2 ) );
               sort01('Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RESULT( 3 ) );
            }

            // Use SBDSQR to form the SVD of the bidiagonal matrix B:
            // B := U * S1 * VT, and compute Z = U' * Y.

            scopy(MNMIN, BD, 1, S1, 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK, 1 );
            slacpy(' ', M, NRHS, Y, LDX, Z, LDX );
            slaset('Full', MNMIN, MNMIN, ZERO, ONE, U, LDPT );
            slaset('Full', MNMIN, MNMIN, ZERO, ONE, VT, LDPT );

            sbdsqr(UPLO, MNMIN, MNMIN, MNMIN, NRHS, S1, WORK, VT, LDPT, U, LDPT, Z, LDX, WORK( MNMIN+1 ), IINFO );

            // Check error code from SBDSQR.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSQR(vects)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[4] = ULPINV;
                  GO TO 270;
               }
            }

            // Use SBDSQR to compute only the singular values of the
            // bidiagonal matrix B;  U, VT, and Z should not be modified.

            scopy(MNMIN, BD, 1, S2, 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK, 1 );

            sbdsqr(UPLO, MNMIN, 0, 0, 0, S2, WORK, VT, LDPT, U, LDPT, Z, LDX, WORK( MNMIN+1 ), IINFO );

            // Check error code from SBDSQR.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSQR(values)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[9] = ULPINV;
                  GO TO 270;
               }
            }

            // Test 4:  Check the decomposition B := U * S1 * VT
                 // 5:  Check the computation Z := U' * Y
                 // 6:  Check the orthogonality of U
                 // 7:  Check the orthogonality of VT

            sbdt03(UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, WORK, RESULT( 4 ) );
            sbdt02(MNMIN, NRHS, Y, LDX, Z, LDX, U, LDPT, WORK, RESULT( 5 ) );
            sort01('Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK, RESULT( 6 ) );
            sort01('Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK, RESULT( 7 ) );

            // Test 8:  Check that the singular values are sorted in
                     // non-increasing order and are non-negative

            RESULT[8] = ZERO;
            for (I = 1; I <= MNMIN - 1; I++) { // 110
               if[S1( I ) < S1( I+1 ) ) RESULT( 8] = ULPINV;
               IF[S1( I ) < ZERO ) RESULT( 8] = ULPINV;
            } // 110
            if ( MNMIN >= 1 ) {
               if[S1( MNMIN ) < ZERO ) RESULT( 8] = ULPINV;
            }

            // Test 9:  Compare SBDSQR with and without singular vectors

            TEMP2 = ZERO;

            for (J = 1; J <= MNMIN; J++) { // 120
               TEMP1 = ABS( S1( J )-S2( J ) ) / max( sqrt( UNFL )*max( S1( 1 ), ONE ), ULP*max( ( S1( J ) ).abs(), ( S2( J ) ) ) ).abs();
               TEMP2 = max( TEMP1, TEMP2 );
            } // 120

            RESULT[9] = TEMP2;

            // Test 10:  Sturm sequence test of singular values
                      // Go up by factors of two until it succeeds

            TEMP1 = THRESH*( HALF-ULP );

            for (J = 0; J <= LOG2UI; J++) { // 130
                // CALL SSVDCH( MNMIN, BD, BE, S1, TEMP1, IINFO )
               if (IINFO == 0) GO TO 140;
               TEMP1 = TEMP1*TWO;
            } // 130

            } // 140
            RESULT[10] = TEMP1;

            // Use SBDSQR to form the decomposition A := (QU) S (VT PT)
            // from the bidiagonal form A := Q B PT.

            if ( !BIDIAG ) {
               scopy(MNMIN, BD, 1, S2, 1 );
               if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK, 1 );

               sbdsqr(UPLO, MNMIN, N, M, NRHS, S2, WORK, PT, LDPT, Q, LDQ, Y, LDX, WORK( MNMIN+1 ), IINFO );

               // Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
                    // 12:  Check the computation Z := U' * Q' * X
                    // 13:  Check the orthogonality of Q*U
                    // 14:  Check the orthogonality of VT*PT

               sbdt01(M, N, 0, A, LDA, Q, LDQ, S2, DUMMA, PT, LDPT, WORK, RESULT( 11 ) );
               sbdt02(M, NRHS, X, LDX, Y, LDX, Q, LDQ, WORK, RESULT( 12 ) );
               sort01('Columns', M, MQ, Q, LDQ, WORK, LWORK, RESULT( 13 ) );
               sort01('Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RESULT( 14 ) );
            }

            // Use SBDSDC to form the SVD of the bidiagonal matrix B:
            // B := U * S1 * VT

            scopy(MNMIN, BD, 1, S1, 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK, 1 );
            slaset('Full', MNMIN, MNMIN, ZERO, ONE, U, LDPT );
            slaset('Full', MNMIN, MNMIN, ZERO, ONE, VT, LDPT );

            sbdsdc(UPLO, 'I', MNMIN, S1, WORK, U, LDPT, VT, LDPT, DUM, IDUM, WORK( MNMIN+1 ), IWORK, IINFO );

            // Check error code from SBDSDC.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSDC(vects)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[15] = ULPINV;
                  GO TO 270;
               }
            }

            // Use SBDSDC to compute only the singular values of the
            // bidiagonal matrix B;  U and VT should not be modified.

            scopy(MNMIN, BD, 1, S2, 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK, 1 );

            sbdsdc(UPLO, 'N', MNMIN, S2, WORK, DUM, 1, DUM, 1, DUM, IDUM, WORK( MNMIN+1 ), IWORK, IINFO );

            // Check error code from SBDSDC.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSDC(values)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[18] = ULPINV;
                  GO TO 270;
               }
            }

            // Test 15:  Check the decomposition B := U * S1 * VT
                 // 16:  Check the orthogonality of U
                 // 17:  Check the orthogonality of VT

            sbdt03(UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, WORK, RESULT( 15 ) );
            sort01('Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK, RESULT( 16 ) );
            sort01('Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK, RESULT( 17 ) );

            // Test 18:  Check that the singular values are sorted in
                      // non-increasing order and are non-negative

            RESULT[18] = ZERO;
            for (I = 1; I <= MNMIN - 1; I++) { // 150
               if[S1( I ) < S1( I+1 ) ) RESULT( 18] = ULPINV;
               IF[S1( I ) < ZERO ) RESULT( 18] = ULPINV;
            } // 150
            if ( MNMIN >= 1 ) {
               if[S1( MNMIN ) < ZERO ) RESULT( 18] = ULPINV;
            }

            // Test 19:  Compare SBDSQR with and without singular vectors

            TEMP2 = ZERO;

            for (J = 1; J <= MNMIN; J++) { // 160
               TEMP1 = ABS( S1( J )-S2( J ) ) / max( sqrt( UNFL )*max( S1( 1 ), ONE ), ULP*max( ( S1( 1 ) ).abs(), ( S2( 1 ) ) ) ).abs();
               TEMP2 = max( TEMP1, TEMP2 );
            } // 160

            RESULT[19] = TEMP2;


            // Use SBDSVDX to compute the SVD of the bidiagonal matrix B:
            // B := U * S1 * VT

            if ( JTYPE == 10 || JTYPE == 16 ) {
               // =================================
               // Matrix types temporarily disabled
               // =================================
               RESULT[20:34] = ZERO;
               GO TO 270;
            }

            IWBS = 1;
            IWBD = IWBS + MNMIN;
            IWBE = IWBD + MNMIN;
            IWBZ = IWBE + MNMIN;
            IWWORK = IWBZ + 2*MNMIN*(MNMIN+1);
            MNMIN2 = max( 1,MNMIN*2 );

            scopy(MNMIN, BD, 1, WORK( IWBD ), 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK( IWBE ), 1 );

            sbdsvdx(UPLO, 'V', 'A', MNMIN, WORK( IWBD ), WORK( IWBE ), ZERO, ZERO, 0, 0, NS1, S1, WORK( IWBZ ), MNMIN2, WORK( IWWORK ), IWORK, IINFO);

            // Check error code from SBDSVDX.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSVDX(vects,A)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[20] = ULPINV;
                  GO TO 270;
               }
            }

            J = IWBZ;
            for (I = 1; I <= NS1; I++) { // 170
               scopy(MNMIN, WORK( J ), 1, U( 1,I ), 1 );
               J = J + MNMIN;
               scopy(MNMIN, WORK( J ), 1, VT( I,1 ), LDPT );
               J = J + MNMIN;
            } // 170

            // Use SBDSVDX to compute only the singular values of the
            // bidiagonal matrix B;  U and VT should not be modified.

            if ( JTYPE == 9 ) {
               // =================================
               // Matrix types temporarily disabled
               // =================================
               RESULT[24] = ZERO;
               GO TO 270;
            }

            scopy(MNMIN, BD, 1, WORK( IWBD ), 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK( IWBE ), 1 );

            sbdsvdx(UPLO, 'N', 'A', MNMIN, WORK( IWBD ), WORK( IWBE ), ZERO, ZERO, 0, 0, NS2, S2, WORK( IWBZ ), MNMIN2, WORK( IWWORK ), IWORK, IINFO );

            // Check error code from SBDSVDX.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSVDX(values,A)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[24] = ULPINV;
                  GO TO 270;
               }
            }

            // Save S1 for tests 30-34.

            scopy(MNMIN, S1, 1, WORK( IWBS ), 1 );

            // Test 20:  Check the decomposition B := U * S1 * VT
                 // 21:  Check the orthogonality of U
                 // 22:  Check the orthogonality of VT
                 // 23:  Check that the singular values are sorted in
                      // non-increasing order and are non-negative
                 // 24:  Compare SBDSVDX with and without singular vectors

            sbdt03(UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, WORK( IWBS+MNMIN ), RESULT( 20 ) );
            sort01('Columns', MNMIN, MNMIN, U, LDPT, WORK( IWBS+MNMIN ), LWORK-MNMIN, RESULT( 21 ) );
            sort01('Rows', MNMIN, MNMIN, VT, LDPT, WORK( IWBS+MNMIN ), LWORK-MNMIN, RESULT( 22) );

            RESULT[23] = ZERO;
            for (I = 1; I <= MNMIN - 1; I++) { // 180
               if[S1( I ) < S1( I+1 ) ) RESULT( 23] = ULPINV;
               IF[S1( I ) < ZERO ) RESULT( 23] = ULPINV;
            } // 180
            if ( MNMIN >= 1 ) {
               if[S1( MNMIN ) < ZERO ) RESULT( 23] = ULPINV;
            }

            TEMP2 = ZERO;
            for (J = 1; J <= MNMIN; J++) { // 190
               TEMP1 = ABS( S1( J )-S2( J ) ) / max( sqrt( UNFL )*max( S1( 1 ), ONE ), ULP*max( ( S1( 1 ) ).abs(), ( S2( 1 ) ) ) ).abs();
               TEMP2 = max( TEMP1, TEMP2 );
            } // 190
            RESULT[24] = TEMP2;
            ANORM = S1( 1 );

            // Use SBDSVDX with RANGE='I': choose random values for IL and
            // IU, and ask for the IL-th through IU-th singular values
            // and corresponding vectors.

            for (I = 1; I <= 4; I++) { // 200
               ISEED2[I] = ISEED( I );
            } // 200
            if ( MNMIN <= 1 ) {
               IL = 1;
               IU = MNMIN;
            } else {
               IL = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) );
               IU = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) );
               if ( IU < IL ) {
                  ITEMP = IU;
                  IU = IL;
                  IL = ITEMP;
               }
            }

            scopy(MNMIN, BD, 1, WORK( IWBD ), 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK( IWBE ), 1 );

            sbdsvdx(UPLO, 'V', 'I', MNMIN, WORK( IWBD ), WORK( IWBE ), ZERO, ZERO, IL, IU, NS1, S1, WORK( IWBZ ), MNMIN2, WORK( IWWORK ), IWORK, IINFO);

            // Check error code from SBDSVDX.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSVDX(vects,I)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[25] = ULPINV;
                  GO TO 270;
               }
            }

            J = IWBZ;
            for (I = 1; I <= NS1; I++) { // 210
               scopy(MNMIN, WORK( J ), 1, U( 1,I ), 1 );
               J = J + MNMIN;
               scopy(MNMIN, WORK( J ), 1, VT( I,1 ), LDPT );
               J = J + MNMIN;
            } // 210

            // Use SBDSVDX to compute only the singular values of the
            // bidiagonal matrix B;  U and VT should not be modified.

            scopy(MNMIN, BD, 1, WORK( IWBD ), 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK( IWBE ), 1 );

            sbdsvdx(UPLO, 'N', 'I', MNMIN, WORK( IWBD ), WORK( IWBE ), ZERO, ZERO, IL, IU, NS2, S2, WORK( IWBZ ), MNMIN2, WORK( IWWORK ), IWORK, IINFO );

            // Check error code from SBDSVDX.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSVDX(values,I)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[29] = ULPINV;
                  GO TO 270;
               }
            }

            // Test 25:  Check S1 - U' * B * VT'
                 // 26:  Check the orthogonality of U
                 // 27:  Check the orthogonality of VT
                 // 28:  Check that the singular values are sorted in
                      // non-increasing order and are non-negative
                 // 29:  Compare SBDSVDX with and without singular vectors

            sbdt04(UPLO, MNMIN, BD, BE, S1, NS1, U, LDPT, VT, LDPT, WORK( IWBS+MNMIN ), RESULT( 25 ) );
            sort01('Columns', MNMIN, NS1, U, LDPT, WORK( IWBS+MNMIN ), LWORK-MNMIN, RESULT( 26 ) );
            sort01('Rows', NS1, MNMIN, VT, LDPT, WORK( IWBS+MNMIN ), LWORK-MNMIN, RESULT( 27 ) );

            RESULT[28] = ZERO;
            for (I = 1; I <= NS1 - 1; I++) { // 220
               if[S1( I ) < S1( I+1 ) ) RESULT( 28] = ULPINV;
               IF[S1( I ) < ZERO ) RESULT( 28] = ULPINV;
            } // 220
            if ( NS1 >= 1 ) {
               if[S1( NS1 ) < ZERO ) RESULT( 28] = ULPINV;
            }

            TEMP2 = ZERO;
            for (J = 1; J <= NS1; J++) { // 230
               TEMP1 = ABS( S1( J )-S2( J ) ) / max( sqrt( UNFL )*max( S1( 1 ), ONE ), ULP*max( ( S1( 1 ) ).abs(), ( S2( 1 ) ) ) ).abs();
               TEMP2 = max( TEMP1, TEMP2 );
            } // 230
            RESULT[29] = TEMP2;

            // Use SBDSVDX with RANGE='V': determine the values VL and VU
            // of the IL-th and IU-th singular values and ask for all
            // singular values in this range.

            scopy(MNMIN, WORK( IWBS ), 1, S1, 1 );

            if ( MNMIN > 0 ) {
               if ( IL != 1 ) {
                  VU = S1( IL ) + max( HALF*ABS( S1( IL )-S1( IL-1 ) ), ULP*ANORM, TWO*RTUNFL );
               } else {
                  VU = S1( 1 ) + max( HALF*ABS( S1( MNMIN )-S1( 1 ) ), ULP*ANORM, TWO*RTUNFL );
               }
               if ( IU != NS1 ) {
                  VL = S1( IU ) - max( ULP*ANORM, TWO*RTUNFL, HALF*ABS( S1( IU+1 )-S1( IU ) ) );
               } else {
                  VL = S1( NS1 ) - max( ULP*ANORM, TWO*RTUNFL, HALF*ABS( S1( MNMIN )-S1( 1 ) ) );
               }
               VL = max( VL,ZERO );
               VU = max( VU,ZERO );
               if (VL >= VU) VU = max( VU*2, VU+VL+HALF );
            } else {
               VL = ZERO;
               VU = ONE;
            }

            scopy(MNMIN, BD, 1, WORK( IWBD ), 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK( IWBE ), 1 );

            sbdsvdx(UPLO, 'V', 'V', MNMIN, WORK( IWBD ), WORK( IWBE ), VL, VU, 0, 0, NS1, S1, WORK( IWBZ ), MNMIN2, WORK( IWWORK ), IWORK, IINFO );

            // Check error code from SBDSVDX.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSVDX(vects,V)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[30] = ULPINV;
                  GO TO 270;
               }
            }

            J = IWBZ;
            for (I = 1; I <= NS1; I++) { // 240
               scopy(MNMIN, WORK( J ), 1, U( 1,I ), 1 );
               J = J + MNMIN;
               scopy(MNMIN, WORK( J ), 1, VT( I,1 ), LDPT );
               J = J + MNMIN;
            } // 240

            // Use SBDSVDX to compute only the singular values of the
            // bidiagonal matrix B;  U and VT should not be modified.

            scopy(MNMIN, BD, 1, WORK( IWBD ), 1 );
            if (MNMIN > 0) scopy( MNMIN-1, BE, 1, WORK( IWBE ), 1 );

            sbdsvdx(UPLO, 'N', 'V', MNMIN, WORK( IWBD ), WORK( IWBE ), VL, VU, 0, 0, NS2, S2, WORK( IWBZ ), MNMIN2, WORK( IWWORK ), IWORK, IINFO );

            // Check error code from SBDSVDX.

            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9998 )'SBDSVDX(values,V)', IINFO, M, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[34] = ULPINV;
                  GO TO 270;
               }
            }

            // Test 30:  Check S1 - U' * B * VT'
                 // 31:  Check the orthogonality of U
                 // 32:  Check the orthogonality of VT
                 // 33:  Check that the singular values are sorted in
                      // non-increasing order and are non-negative
                 // 34:  Compare SBDSVDX with and without singular vectors

            sbdt04(UPLO, MNMIN, BD, BE, S1, NS1, U, LDPT, VT, LDPT, WORK( IWBS+MNMIN ), RESULT( 30 ) );
            sort01('Columns', MNMIN, NS1, U, LDPT, WORK( IWBS+MNMIN ), LWORK-MNMIN, RESULT( 31 ) );
            sort01('Rows', NS1, MNMIN, VT, LDPT, WORK( IWBS+MNMIN ), LWORK-MNMIN, RESULT( 32 ) );

            RESULT[33] = ZERO;
            for (I = 1; I <= NS1 - 1; I++) { // 250
               if[S1( I ) < S1( I+1 ) ) RESULT( 28] = ULPINV;
               IF[S1( I ) < ZERO ) RESULT( 28] = ULPINV;
            } // 250
            if ( NS1 >= 1 ) {
               if[S1( NS1 ) < ZERO ) RESULT( 28] = ULPINV;
            }

            TEMP2 = ZERO;
            for (J = 1; J <= NS1; J++) { // 260
               TEMP1 = ABS( S1( J )-S2( J ) ) / max( sqrt( UNFL )*max( S1( 1 ), ONE ), ULP*max( ( S1( 1 ) ).abs(), ( S2( 1 ) ) ) ).abs();
               TEMP2 = max( TEMP1, TEMP2 );
            } // 260
            RESULT[34] = TEMP2;

            // End of Loop -- Check for RESULT(j) > THRESH

            } // 270

            for (J = 1; J <= 34; J++) { // 280
               if ( RESULT( J ) >= THRESH ) {
                  if (NFAIL == 0) slahd2( NOUT, PATH );
                  WRITE( NOUT, FMT = 9999 )M, N, JTYPE, IOLDSD, J, RESULT( J );
                  NFAIL = NFAIL + 1;
               }
            } // 280
            if ( !BIDIAG ) {
               NTEST = NTEST + 34;
            } else {
               NTEST = NTEST + 30;
            }

         } // 290
      } // 300

      // Summary

      alasum(PATH, NOUT, NFAIL, NTEST, 0 );

      return;
 9999 FORMAT( ' M=', I5, ', N=', I5, ', type ', I2, ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 );
 9998 FORMAT( ' SCHKBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

      }
