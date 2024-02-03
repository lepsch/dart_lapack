      SUBROUTINE SLASRT( ID, N, D, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             ID;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                SELECT;
      const              SELECT = 20 ;
      // ..
      // .. Local Scalars ..
      int                DIR, ENDD, I, J, START, STKPNT;
      REAL               D1, D2, D3, DMNMX, TMP
      // ..
      // .. Local Arrays ..
      int                STACK( 2, 32 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      DIR = -1
      if ( LSAME( ID, 'D' ) ) {
         DIR = 0
      } else if ( LSAME( ID, 'I' ) ) {
         DIR = 1
      }
      if ( DIR.EQ.-1 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SLASRT', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.LE.1 ) RETURN

      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      if ( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) {

         // Do Insertion sort on D( START:ENDD )

         if ( DIR.EQ.0 ) {

            // Sort into decreasing order

            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  if ( D( J ).GT.D( J-1 ) ) {
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  } else {
                     GO TO 30
                  }
   20          CONTINUE
   30       CONTINUE

         } else {

            // Sort into increasing order

            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  if ( D( J ).LT.D( J-1 ) ) {
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  } else {
                     GO TO 50
                  }
   40          CONTINUE
   50       CONTINUE

         }

      } else if ( ENDD-START.GT.SELECT ) {

         // Partition D( START:ENDD ) and stack parts, largest one first

         // Choose partition entry as median of 3

         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         if ( D1.LT.D2 ) {
            if ( D3.LT.D1 ) {
               DMNMX = D1
            } else if ( D3.LT.D2 ) {
               DMNMX = D3
            } else {
               DMNMX = D2
            }
         } else {
            if ( D3.LT.D2 ) {
               DMNMX = D2
            } else if ( D3.LT.D1 ) {
               DMNMX = D3
            } else {
               DMNMX = D1
            }
         }

         if ( DIR.EQ.0 ) {

            // Sort into decreasing order

            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX ) GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX ) GO TO 80
            if ( I.LT.J ) {
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            }
            if ( J-START.GT.ENDD-J-1 ) {
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            } else {
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            }
         } else {

            // Sort into increasing order

            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX ) GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX ) GO TO 110
            if ( I.LT.J ) {
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            }
            if ( J-START.GT.ENDD-J-1 ) {
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            } else {
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            }
         }
      }
      IF( STKPNT.GT.0 ) GO TO 10
      RETURN

      // End of SLASRT

      }
