Das Programm quad_opt wurde entwickelt um mittels Quadratischer Optimierung aus Messdaten die in Koordinatenform vorliegen, einen passenden Graphen zu errechnen.

Quad_op wurde f�r Linux entwickelt und kann teilweise sowohl direkt via Temrminal als auch mit einem rudiment�ren Men� bedient werden.

Die Messdaten m�ssen in einem txt-File vorliegen deren Daten in folgender Formatierung gespeichert sein m�ssen. In der ersten Zeile muss die Anzahl der Punkte stehen, danach folgen jeweils die x und y Werte
z.B. messung1.txt:
9

0.5 0.5

1.0 1.5

1.5 2.0

2.5 3.0

3.0 3.0

3.5 2.5

4.0 2.0

4.5 1.0

5.0 0.0

Der direktee Aufruf aus dem Terminal 
	Der Aufruf erfolgt mit dem Namen des Messdatenfiles und mit dem Grad des Polynoms das errechnet werden soll: ./quad_opt <filename>.txt <2-9>
	Es k�nnen Polynome mit dem Grad 2 bis 9 genutzt werden.
	Es wird das File <filename>_result.txt erstellt mit den Koeffizienten des Polynoms. 

Wird quad_opt gestartet kann man aus 4 Punkten ausw�hlen:
	(1) "generate testfile" - es soll ein File erstellt werden das Messpunkte enth�lt die ann�hernd nach einem Polynom verteilt sind
	hier wird zun�chst der Name des zu erstellen Testfiles, dann der Grad des Polynoms und die Anzahl der Punkte abgefragt
	anschlie�end wird das File <filename>.txt erstellt

	(2) "load file and opotimize!" - errechnet die Koeffizienten aus den �bergebenen File und Grad des Polynoms
	Es wird das File <filename>_result.txt erstellt das die errechneten Koeffizienten des Polynoms enth�lt. 

	(3) "GnuPlot" - aus dem Messdaten werden die Koeffizienten errechnet, anschlie�end wird das Gnu Scriptfile <filename>.gnu erstellt und aufgerufen
	Das Scriptfile erstellt die Bilddatei <filename>.png in dem die Messpunkte sowie der errechnete Graph dargestellt werden.

	(4) "create 'btest' ... " ruft alle vorhergehenden Punkte nacheinander auf, als Filname ist btest und der Grad 2 eingestellt.

