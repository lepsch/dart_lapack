      SUBROUTINE DDRVRF2( NOUT, NN, NVAL, A, LDA, ARF, AP, ASAV  )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, NN, NOUT;
      // ..
      // .. Array Arguments ..
      int                NVAL( NN );
      double             A( LDA, * ), ARF( * ), AP(*), ASAV( LDA, * );
      // ..

*  =====================================================================
      // ..
      // .. Local Scalars ..
      bool               LOWER, OK1, OK2;
      String             UPLO, CFORM;
      int                I, IFORM, IIN, INFO, IUPLO, J, N, NERRS, NRUN;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      // ..
      // .. External Functions ..
      double             DLARND;
      // EXTERNAL DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTFTTR, DTFTTP, DTRTTF, DTRTTP, DTPTTR, DTPTTF
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      DATA               FORMS / 'N', 'T' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      NRUN = 0
      NERRS = 0
      INFO = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      DO 120 IIN = 1, NN

         N = NVAL( IIN )

         // Do first for UPLO = 'U', then for UPLO = 'L'

         DO 110 IUPLO = 1, 2

            UPLO = UPLOS( IUPLO )
            LOWER = .TRUE.
            IF ( IUPLO.EQ.1 ) LOWER = .FALSE.

            // Do first for CFORM = 'N', then for CFORM = 'T'

            DO 100 IFORM = 1, 2

               CFORM = FORMS( IFORM )

               NRUN = NRUN + 1

               DO J = 1, N
                  DO I = 1, N
                     A( I, J) = DLARND( 2, ISEED )
                  END DO
               END DO

               SRNAMT = 'DTRTTF'
               CALL DTRTTF( CFORM, UPLO, N, A, LDA, ARF, INFO )

               SRNAMT = 'DTFTTP'
               CALL DTFTTP( CFORM, UPLO, N, ARF, AP, INFO )

               SRNAMT = 'DTPTTR'
               CALL DTPTTR( UPLO, N, AP, ASAV, LDA, INFO )

               OK1 = .TRUE.
               if ( LOWER ) {
                  DO J = 1, N
                     DO I = J, N
                        if ( A(I,J).NE.ASAV(I,J) ) {
                           OK1 = .FALSE.
                        }
                     END DO
                  END DO
               } else {
                  DO J = 1, N
                     DO I = 1, J
                        if ( A(I,J).NE.ASAV(I,J) ) {
                           OK1 = .FALSE.
                        }
                     END DO
                  END DO
               }

               NRUN = NRUN + 1

               SRNAMT = 'DTRTTP'
               CALL DTRTTP( UPLO, N, A, LDA, AP, INFO )

               SRNAMT = 'DTPTTF'
               CALL DTPTTF( CFORM, UPLO, N, AP, ARF, INFO )

               SRNAMT = 'DTFTTR'
               CALL DTFTTR( CFORM, UPLO, N, ARF, ASAV, LDA, INFO )

               OK2 = .TRUE.
               if ( LOWER ) {
                  DO J = 1, N
                     DO I = J, N
                        if ( A(I,J).NE.ASAV(I,J) ) {
                           OK2 = .FALSE.
                        }
                     END DO
                  END DO
               } else {
                  DO J = 1, N
                     DO I = 1, J
                        if ( A(I,J).NE.ASAV(I,J) ) {
                           OK2 = .FALSE.
                        }
                     END DO
                  END DO
               }

               if (( .NOT.OK1 ).OR.( .NOT.OK2 )) {
                  if ( NERRS.EQ.0 ) {
                     WRITE( NOUT, * )
                     WRITE( NOUT, FMT = 9999 )
                  }
                  WRITE( NOUT, FMT = 9998 ) N, UPLO, CFORM
                  NERRS = NERRS + 1
               }

  100       CONTINUE
  110    CONTINUE
  120 CONTINUE

      // Print a summary of the results.

      if ( NERRS.EQ.0 ) {
         WRITE( NOUT, FMT = 9997 ) NRUN
      } else {
         WRITE( NOUT, FMT = 9996 ) NERRS, NRUN
      }

 9999 FORMAT( 1X, ' *** Error(s) while testing the RFP conversion',
     +         ' routines ***')
 9998 FORMAT( 1X, '     Error in RFP,conversion routines N=',I5,
     +        ' UPLO=''', A1, ''', FORM =''',A1,'''')
 9997 FORMAT( 1X, 'All tests for the RFP conversion routines passed ( ',
     +        I5,' tests run)')
 9996 FORMAT( 1X, 'RFP conversion routines: ',I5,' out of ',I5,
     +        ' error message recorded')

      RETURN

      // End of DDRVRF2

      }
