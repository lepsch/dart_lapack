      SUBROUTINE ZCHKBD( NSIZES, MVAL, NVAL, NTYPES, DOTYPE, NRHS, ISEED, THRESH, A, LDA, BD, BE, S1, S2, X, LDX, Y, Z, Q, LDQ, PT, LDPT, U, VT, WORK, LWORK, RWORK, NOUT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDPT, LDQ, LDX, LWORK, NOUT, NRHS, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), MVAL( * ), NVAL( * );
      double             BD( * ), BE( * ), RWORK( * ), S1( * ), S2( * );
      COMPLEX*16         A( LDA, * ), PT( LDPT, * ), Q( LDQ, * ), U( LDPT, * ), VT( LDPT, * ), WORK( * ), X( LDX, * ), Y( LDX, * ), Z( LDX, * )
      // ..

* ======================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, HALF;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      int                MAXTYP;
      const              MAXTYP = 16 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN, BIDIAG;
      String             UPLO;
      String             PATH;
      int                I, IINFO, IMODE, ITYPE, J, JCOL, JSIZE, JTYPE, LOG2UI, M, MINWRK, MMAX, MNMAX, MNMIN, MQ, MTYPES, N, NFAIL, NMAX, NTEST;
      double             AMNINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
      int                IOLDSD( 4 ), IWORK( 1 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      double             DUMMA( 1 ), RESULT( 14 );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLARND;
      // EXTERNAL DLAMCH, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASUM, DCOPY, DLAHD2, DSVDCH, XERBLA, ZBDSQR, ZBDT01, ZBDT02, ZBDT03, ZGEBRD, ZGEMM, ZLACPY, ZLASET, ZLATMR, ZLATMS, ZUNGBR, ZUNT01
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
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA            KTYPE / 1, 2, 5*4, 5*6, 3*9, 10 /
      DATA            KMAGN / 2*1, 3*1, 2, 3, 3*1, 2, 3, 1, 2, 3, 0 /
      DATA            KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 0 /
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0

      BADMM = .FALSE.
      BADNN = .FALSE.
      MMAX = 1
      NMAX = 1
      MNMAX = 1
      MINWRK = 1
      DO 10 J = 1, NSIZES
         MMAX = MAX( MMAX, MVAL( J ) )
         IF( MVAL( J ).LT.0 ) BADMM = .TRUE.
         NMAX = MAX( NMAX, NVAL( J ) )
         IF( NVAL( J ).LT.0 ) BADNN = .TRUE.
         MNMAX = MAX( MNMAX, MIN( MVAL( J ), NVAL( J ) ) )
         MINWRK = MAX( MINWRK, 3*( MVAL( J )+NVAL( J ) ), MVAL( J )*( MVAL( J )+MAX( MVAL( J ), NVAL( J ), NRHS )+1 )+NVAL( J )*MIN( NVAL( J ), MVAL( J ) ) )
   10 CONTINUE

      // Check for errors

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADMM ) {
         INFO = -2
      } else if ( BADNN ) {
         INFO = -3
      } else if ( NTYPES.LT.0 ) {
         INFO = -4
      } else if ( NRHS.LT.0 ) {
         INFO = -6
      } else if ( LDA.LT.MMAX ) {
         INFO = -11
      } else if ( LDX.LT.MMAX ) {
         INFO = -17
      } else if ( LDQ.LT.MMAX ) {
         INFO = -21
      } else if ( LDPT.LT.MNMAX ) {
         INFO = -23
      } else if ( MINWRK.GT.LWORK ) {
         INFO = -27
      }

      if ( INFO.NE.0 ) {
         xerbla('ZCHKBD', -INFO );
         RETURN
      }

      // Initialize constants

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'BD'
      NFAIL = 0
      NTEST = 0
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      LOG2UI = INT( LOG( ULPINV ) / LOG( TWO ) )
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
      INFOT = 0

      // Loop over sizes, types

      DO 180 JSIZE = 1, NSIZES
         M = MVAL( JSIZE )
         N = NVAL( JSIZE )
         MNMIN = MIN( M, N )
         AMNINV = ONE / MAX( M, N, 1 )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         DO 170 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 170

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            DO 30 J = 1, 14
               RESULT( J ) = -ONE
   30       CONTINUE

            UPLO = ' '

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

            IF( MTYPES.GT.MAXTYP ) GO TO 100

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE )

   40       CONTINUE
            ANORM = ONE
            GO TO 70

   50       CONTINUE
            ANORM = ( RTOVFL*ULP )*AMNINV
            GO TO 70

   60       CONTINUE
            ANORM = RTUNFL*MAX( M, N )*ULPINV
            GO TO 70

   70       CONTINUE

            zlaset('Full', LDA, N, CZERO, CZERO, A, LDA );
            IINFO = 0
            COND = ULPINV

            BIDIAG = .FALSE.
            if ( ITYPE.EQ.1 ) {

               // Zero matrix

               IINFO = 0

            } else if ( ITYPE.EQ.2 ) {

               // Identity

               DO 80 JCOL = 1, MNMIN
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE

            } else if ( ITYPE.EQ.4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               zlatms(MNMIN, MNMIN, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE.EQ.5 ) {

               // Symmetric, eigenvalues specified

               zlatms(MNMIN, MNMIN, 'S', ISEED, 'S', RWORK, IMODE, COND, ANORM, M, N, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE.EQ.6 ) {

               // Nonsymmetric, singular values specified

               zlatms(M, N, 'S', ISEED, 'N', RWORK, IMODE, COND, ANORM, M, N, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE.EQ.7 ) {

               // Diagonal, random entries

               zlatmr(MNMIN, MNMIN, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( 2*MNMIN+1 ), 1, ONE, 'N', IWORK, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE.EQ.8 ) {

               // Symmetric, random entries

               zlatmr(MNMIN, MNMIN, 'S', ISEED, 'S', WORK, 6, ONE, CONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( M+MNMIN+1 ), 1, ONE, 'N', IWORK, M, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE.EQ.9 ) {

               // Nonsymmetric, random entries

               zlatmr(M, N, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( M+MNMIN+1 ), 1, ONE, 'N', IWORK, M, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE.EQ.10 ) {

               // Bidiagonal, random entries

               TEMP1 = -TWO*LOG( ULP )
               DO 90 J = 1, MNMIN
                  BD( J ) = EXP( TEMP1*DLARND( 2, ISEED ) )
                  IF( J.LT.MNMIN ) BE( J ) = EXP( TEMP1*DLARND( 2, ISEED ) )
   90          CONTINUE

               IINFO = 0
               BIDIAG = .TRUE.
               if ( M.GE.N ) {
                  UPLO = 'U'
               } else {
                  UPLO = 'L'
               }
            } else {
               IINFO = 1
            }

            if ( IINFO.EQ.0 ) {

               // Generate Right-Hand Side

               if ( BIDIAG ) {
                  zlatmr(MNMIN, NRHS, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( MNMIN+1 ), 1, ONE, WORK( 2*MNMIN+1 ), 1, ONE, 'N', IWORK, MNMIN, NRHS, ZERO, ONE, 'NO', Y, LDX, IWORK, IINFO );
               } else {
                  zlatmr(M, NRHS, 'S', ISEED, 'N', WORK, 6, ONE, CONE, 'T', 'N', WORK( M+1 ), 1, ONE, WORK( 2*M+1 ), 1, ONE, 'N', IWORK, M, NRHS, ZERO, ONE, 'NO', X, LDX, IWORK, IINFO );
               }
            }

            // Error Exit

            if ( IINFO.NE.0 ) {
               WRITE( NOUT, FMT = 9998 )'Generator', IINFO, M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

  100       CONTINUE

            // Call ZGEBRD and ZUNGBR to compute B, Q, and P, do tests.

            if ( .NOT.BIDIAG ) {

               // Compute transformations to reduce A to bidiagonal form:
               // B := Q' * A * P.

               zlacpy(' ', M, N, A, LDA, Q, LDQ );
               zgebrd(M, N, Q, LDQ, BD, BE, WORK, WORK( MNMIN+1 ), WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO );

               // Check error code from ZGEBRD.

               if ( IINFO.NE.0 ) {
                  WRITE( NOUT, FMT = 9998 )'ZGEBRD', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               zlacpy(' ', M, N, Q, LDQ, PT, LDPT );
               if ( M.GE.N ) {
                  UPLO = 'U'
               } else {
                  UPLO = 'L'
               }

               // Generate Q

               MQ = M
               IF( NRHS.LE.0 ) MQ = MNMIN                CALL ZUNGBR( 'Q', M, MQ, N, Q, LDQ, WORK, WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO )

               // Check error code from ZUNGBR.

               if ( IINFO.NE.0 ) {
                  WRITE( NOUT, FMT = 9998 )'ZUNGBR(Q)', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Generate P'

               zungbr('P', MNMIN, N, M, PT, LDPT, WORK( MNMIN+1 ), WORK( 2*MNMIN+1 ), LWORK-2*MNMIN, IINFO );

               // Check error code from ZUNGBR.

               if ( IINFO.NE.0 ) {
                  WRITE( NOUT, FMT = 9998 )'ZUNGBR(P)', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Apply Q' to an M by NRHS matrix X:  Y := Q' * X.

               zgemm('Conjugate transpose', 'No transpose', M, NRHS, M, CONE, Q, LDQ, X, LDX, CZERO, Y, LDX );

               // Test 1:  Check the decomposition A := Q * B * PT
                    // 2:  Check the orthogonality of Q
                    // 3:  Check the orthogonality of PT

               zbdt01(M, N, 1, A, LDA, Q, LDQ, BD, BE, PT, LDPT, WORK, RWORK, RESULT( 1 ) )                CALL ZUNT01( 'Columns', M, MQ, Q, LDQ, WORK, LWORK, RWORK, RESULT( 2 ) )                CALL ZUNT01( 'Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RWORK, RESULT( 3 ) );
            }

            // Use ZBDSQR to form the SVD of the bidiagonal matrix B:
            // B := U * S1 * VT, and compute Z = U' * Y.

            dcopy(MNMIN, BD, 1, S1, 1 );
            IF( MNMIN.GT.0 ) CALL DCOPY( MNMIN-1, BE, 1, RWORK, 1 )
            zlacpy(' ', M, NRHS, Y, LDX, Z, LDX );
            zlaset('Full', MNMIN, MNMIN, CZERO, CONE, U, LDPT );
            zlaset('Full', MNMIN, MNMIN, CZERO, CONE, VT, LDPT );

            zbdsqr(UPLO, MNMIN, MNMIN, MNMIN, NRHS, S1, RWORK, VT, LDPT, U, LDPT, Z, LDX, RWORK( MNMIN+1 ), IINFO );

            // Check error code from ZBDSQR.

            if ( IINFO.NE.0 ) {
               WRITE( NOUT, FMT = 9998 )'ZBDSQR(vects)', IINFO, M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 4 ) = ULPINV
                  GO TO 150
               }
            }

            // Use ZBDSQR to compute only the singular values of the
            // bidiagonal matrix B;  U, VT, and Z should not be modified.

            dcopy(MNMIN, BD, 1, S2, 1 );
            IF( MNMIN.GT.0 ) CALL DCOPY( MNMIN-1, BE, 1, RWORK, 1 )

            zbdsqr(UPLO, MNMIN, 0, 0, 0, S2, RWORK, VT, LDPT, U, LDPT, Z, LDX, RWORK( MNMIN+1 ), IINFO );

            // Check error code from ZBDSQR.

            if ( IINFO.NE.0 ) {
               WRITE( NOUT, FMT = 9998 )'ZBDSQR(values)', IINFO, M, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 9 ) = ULPINV
                  GO TO 150
               }
            }

            // Test 4:  Check the decomposition B := U * S1 * VT
                 // 5:  Check the computation Z := U' * Y
                 // 6:  Check the orthogonality of U
                 // 7:  Check the orthogonality of VT

            zbdt03(UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, WORK, RESULT( 4 ) )             CALL ZBDT02( MNMIN, NRHS, Y, LDX, Z, LDX, U, LDPT, WORK, RWORK, RESULT( 5 ) )             CALL ZUNT01( 'Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK, RWORK, RESULT( 6 ) )             CALL ZUNT01( 'Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK, RWORK, RESULT( 7 ) );

            // Test 8:  Check that the singular values are sorted in
                     // non-increasing order and are non-negative

            RESULT( 8 ) = ZERO
            DO 110 I = 1, MNMIN - 1
               IF( S1( I ).LT.S1( I+1 ) ) RESULT( 8 ) = ULPINV                IF( S1( I ).LT.ZERO ) RESULT( 8 ) = ULPINV
  110       CONTINUE
            if ( MNMIN.GE.1 ) {
               IF( S1( MNMIN ).LT.ZERO ) RESULT( 8 ) = ULPINV
            }

            // Test 9:  Compare ZBDSQR with and without singular vectors

            TEMP2 = ZERO

            DO 120 J = 1, MNMIN
               TEMP1 = ABS( S1( J )-S2( J ) ) / MAX( SQRT( UNFL )*MAX( S1( 1 ), ONE ), ULP*MAX( ABS( S1( J ) ), ABS( S2( J ) ) ) )
               TEMP2 = MAX( TEMP1, TEMP2 )
  120       CONTINUE

            RESULT( 9 ) = TEMP2

            // Test 10:  Sturm sequence test of singular values
                      // Go up by factors of two until it succeeds

            TEMP1 = THRESH*( HALF-ULP )

            DO 130 J = 0, LOG2UI
               dsvdch(MNMIN, BD, BE, S1, TEMP1, IINFO );
               IF( IINFO.EQ.0 ) GO TO 140
               TEMP1 = TEMP1*TWO
  130       CONTINUE

  140       CONTINUE
            RESULT( 10 ) = TEMP1

            // Use ZBDSQR to form the decomposition A := (QU) S (VT PT)
            // from the bidiagonal form A := Q B PT.

            if ( .NOT.BIDIAG ) {
               dcopy(MNMIN, BD, 1, S2, 1 );
               IF( MNMIN.GT.0 ) CALL DCOPY( MNMIN-1, BE, 1, RWORK, 1 )

               zbdsqr(UPLO, MNMIN, N, M, NRHS, S2, RWORK, PT, LDPT, Q, LDQ, Y, LDX, RWORK( MNMIN+1 ), IINFO );

               // Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
                    // 12:  Check the computation Z := U' * Q' * X
                    // 13:  Check the orthogonality of Q*U
                    // 14:  Check the orthogonality of VT*PT

               zbdt01(M, N, 0, A, LDA, Q, LDQ, S2, DUMMA, PT, LDPT, WORK, RWORK, RESULT( 11 ) )                CALL ZBDT02( M, NRHS, X, LDX, Y, LDX, Q, LDQ, WORK, RWORK, RESULT( 12 ) )                CALL ZUNT01( 'Columns', M, MQ, Q, LDQ, WORK, LWORK, RWORK, RESULT( 13 ) )                CALL ZUNT01( 'Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RWORK, RESULT( 14 ) );
            }

            // End of Loop -- Check for RESULT(j) > THRESH

  150       CONTINUE
            DO 160 J = 1, 14
               if ( RESULT( J ).GE.THRESH ) {
                  IF( NFAIL.EQ.0 ) CALL DLAHD2( NOUT, PATH )                   WRITE( NOUT, FMT = 9999 )M, N, JTYPE, IOLDSD, J, RESULT( J )
                  NFAIL = NFAIL + 1
               }
  160       CONTINUE
            if ( .NOT.BIDIAG ) {
               NTEST = NTEST + 14
            } else {
               NTEST = NTEST + 5
            }

  170    CONTINUE
  180 CONTINUE

      // Summary

      alasum(PATH, NOUT, NFAIL, NTEST, 0 );

      RETURN

      // End of ZCHKBD

 9999 FORMAT( ' M=', I5, ', N=', I5, ', type ', I2, ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9998 FORMAT( ' ZCHKBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      }
