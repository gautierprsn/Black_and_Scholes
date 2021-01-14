#pragma once
#ifndef Resolv_H
#define ResolvE_H

/*!
 * \file methode_asbtrait.h
 * \brief Creation des differentes méthodes de résolution
 * \author Gautier Poursin,Théo Lazzaroni
 */


 /*! \namespace ensiie
  *
  * esage de nommage regroupant tous le projet
  */

namespace ensiie {
	class Resolv 
	{
	protected:
		float dT; /*!< discrétisation du temps*/
		float dL;/*!< discrétisation de l'espace*/

	public:
		/*!
	 *  \brief Destructeur
	 *
	 *  Destructeur de la classe Resolv
	 */
		~Resolv() {};
		/*!
	 *  \brief Constructeur
	 *
	 *  Constructeur de la classe Resolv
	 *
	 *  \param float: T le temps
	 *  \param float: L l'espace
	 *  \param int: N nombre d'interval pour l'espace
	 *   \param int M: nombre d'interval pour le temps
	 */
		Resolv(float, float,int, int);

		/*!
	 *  \brief création des coefficients de la matrice
	 *
	 *  Methode qui permet de calculer les coefficients non nuls des matrices
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice
	 * \param float: sigma donné dans l'énoncé
	 * \param float: r donné dans l'énoncé

	 */
		virtual void computeCoef(float*, float*, float*, float*,int,float, float);

		/*!
	 *  \brief création des matrices
	 *
	 *  Permet de créer les matrices avec les coefficients correspondants et calculés précédement
	 * 
	 *  \param float**: matrice vide
	 * \param float**: matrice vide
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice
	 */
		virtual void computeMatrix(float** , float**,float*,float*,float*,float*,int);

		/*!
	 *  \brief fonction d'affichage du temps et de l'espace
	 *
	 * \param int: M nombre d'interval
	 * \param int: N nombre d'interval
	 * \param float*: V vecteur calculé à chaque étape
	 * 
	 *  Affiche le vecteur calculé
	 */
		virtual void print(int,int,float*) const=0;

		/*!
	 *  \brief fonction permettant d'avoir dT correspondant
	 * \return float 
	 */
		virtual float get_dT() { return dT; }

		/*!
	 *  \brief fonction permettant d'avoir dL correspondant
	 * \return float
	 */
		virtual float get_dL() { return dL; }

		/*!
	 *  \brief fonction permettant d'avoir l'id de la méthode: 1 pour un call, 0 pour un put
	 * \return int
	 */
		virtual int get_id() const = 0;

		/*!
	 *  \brief fonction permettant d'initialiser le vecteur V
	 * \param int: id de la méthode
	 * \param float*: V le vecteur à initialiser
	 *  \param int: M la taille du vecteur
	 * \param float: K la constante de l'énoncé
	 * \return float*
	 */
		virtual float* initialize(int,float*,int, float) const = 0;

		/*!
	 *  \brief fonction permettant de calculer la puissance d'un  nombre
	 * \param float : a 
	 * \param int: b la puissance
	 * \return float
	 */
		virtual float pow(float a, int b) {
			if (b == 0) return 1;
			else return a * pow(a, b - 1);
		}

		/*!
	 *  \brief fonction permettant d'avoir le maximum entre 2 nombre
	 * \param float : a
	 * \param float: b
	 * \return float
	 */
		virtual float max(float a, float b) {
			if (a >= b) return a;
			return b;
		}

		/*!
	 *  \brief fonction permettant d'avoir le minimum entre 2 nombre
	 * \param float : a
	 * \param float: b
	 * \return float
	 */
		virtual float min(float a, float b) const {
			if (a <= b) return a;
			return b;
		}


		/*!
	 *  \brief libération de la mémoire prise
	 *
	 *  Methode qui permet de libérer la mémoire utilisée
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param float**: matrice M1
	 * \param float**: matrice M2
	 * \param int: M
	 * \param float*: Vecteur V
	 * \param float*: vecteur k

	 */
		virtual void deleted(float* a, float* b, float* c, float* d, float** M1, float** M2,int M,float*V,float*k) {
			delete[]a;
			delete[]b;
			delete[]c;
			delete[]d;

			for (int i = 0; i < M - 1; i++) {
				delete M1[i];
				delete M2[i];
			}

			delete[]M1;
			delete[]M2;
			delete[]V;
			delete[]k;

		}    

		/*!
	 *  \brief méthode déocmpoisition LU
	 *
	 *  Methode qui permet de mettre en place la décomposition LU pour rédoure Ax=B avec A tridiagonale
	 *
	 * \param float**: matrice A
	 * \param float*: vecteur B
	 * \param int: n

	 */
		virtual float * decompLU(float**& , float*, int);

		/*!
	 *  \brief méthode produit matrice par vecteur
	 *
	 *  Methode qui permet de calculer A*x avec A matrice et x vecteur
	 *
	 * \param float**: matrice A
	 * \param float*: vecteur x
	 * \param int : n nombre de ligne de la matrice
	 * \param int: m nombre de colonne de la matrice

	 */
	virtual float* product(float**&, float*, int,int);

	/*!
	 *  \brief met a jour le vecteur V a chaque étape
	 *
	 *  Methode qui permet de calculer le vecteur Vn en fonction de Vn+1
	 *
	 * \param float**: matrice M1
	 * \param float**: matrice M2
	 * \param float*: vecteur V
	 * \param int: M nombre de lignes de M1 et M2
	 * \param int: N nombre de colonnes de M1 et M2
	 * \param float*: coefficients de la ligne sous la diagonale de M1
	 * \param float*: vecteur k
	 * \param float*: vector1 vecteur temporaire
	 * \param float: K donné
	 * \param float: r donné
	 * \param float: T donné

	 */
	virtual void calculate(float**&  , float**&  , float*, int, int, float*, float*, float*, float, float, float) =0;

		//SCHEMA DIFF FINIES

		/*!
	 *  \brief création des matrices pour la methode des différences finies
	 *
	 *  Permet de créer la matrice avec les coefficients correspondants 
	 *
	 *  \param float**: matrice nulle
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice carrée
	 */
		virtual void computeMatrixDiffFinies(float**  , float*, float*, int );


		/*!
	 *  \brief création des coefficients de la matrice pour les differences finies
	 *
	 *  Methode qui permet de calculer les coefficients non nuls des matrices de la méthode des différences finies
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice carrée
	 * \param float: mu prix égal a 1
	 * \param float: sigma

	 */
		virtual void computeCoefDiffFinies(float*, float*, int, float,float);

		/*!
	 *  \brief fonction permettant d'avoir dT correspondant pour la méthode des différences finies
	 * \return double
	 */
		virtual double get_dT_tilde(float sigma) { return dT*0.5*sigma*sigma; }

		/*!
	 *  \brief fonction permettant d'initialiser le vecteur V
	 * \param int: id de la méthode
	 * \param float*: V le vecteur à initialiser
	 *  \param int: M la taille du vecteur
	 * \param float: K la constante de l'énoncé
	 * \return float*
	 */
		virtual float* initializeDiffFinies(int, float*, int, float) const = 0;

		/*!
	 *  \brief met a jour le vecteur V a chaque étape
	 *
	 *  Methode qui permet de calculer le vecteur Vn en fonction de Vn+1
	 *
	 * \param float**: matrice M1
	 * \param float*: vecteur V
	 * \param int: M nombre de lignes de M1 et M2
	 * \param int: N nombre de colonnes de M1 et M2
	 * \param float*: coefficients de la ligne sous la diagonale de M1
	 * \param float*: vecteur k
	 * \param float*: vector1 vecteur temporaire
	 * \param float: K donné
	 * \param float: r donné
	 * \param float: T donné
	 * \param float: sigma donné

	 */
		virtual void calculateDiffFinies(float**& , float*, int, int, float*, float*, float*, float, float, float, float)=0;

		/*!
	 *  \brief libération de la mémoire prise
	 *
	 *  Methode qui permet de libérer la mémoire utilisée
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param float**: matrice M1
	 * \param float**: matrice M2
	 * \param int: M
	 * \param float*: Vecteur V
	 * \param float*: vecteur k

	 */
		virtual void deletedDiffFinies(float* a, float* b,  float** M1, float** M2, int M, float* V, float* k) {
			delete[]a;
			delete[]b;

			for (int i = 0; i < M - 1; i++) {
				delete M1[i];
				delete M2[i];
			}

			delete[]M1;
			delete[]M2;
			delete[]V;
			delete[]k;

		}

	};

	class Put: public ensiie::Resolv {
	private:
		int id;/*!< id de la méthode*/
	public: 
		/*!
	 *  \brief Constructeur
	 *
	 *  Constructeur de la classe Put
	 *
	 *  \param float: T le temps
	 *  \param float: L l'espace
	 *  \param int: N nombre d'interval pour l'espace
	 *   \param int M: nombre d'interval pour le temps
	 * \param int: id de la méthode
	 */
		Put(float, float, int, int, int );

		/*!
	 *  \brief Destructeur
	 *
	 *  Destructeur de la classe Put
	 */
		~Put() {};


		/*!
	 *  \brief création des coefficients de la matrice
	 *
	 *  Methode qui permet de calculer les coefficients non nuls des matrices
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice
	 * \param float: sigma donné dans l'énoncé
	 * \param float: r donné dans l'énoncé

	 */

		//Schema Cranck Nicholsosn
		void computeCoef(float*, float*, float*, float*, int, float, float);

		/*!
	 *  \brief création des matrices
	 *
	 *  Permet de créer les matrices avec les coefficients correspondants et calculés précédement
	 *
	 *  \param float**: matrice vide
	 * \param float**: matrice vide
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice
	 */
		void computeMatrix(float**, float**, float*, float*, float*, float*, int);

		/*!
	 *  \brief fonction permettant d'initialiser le vecteur V
	 * \param int: id de la méthode
	 * \param float*: V le vecteur à initialiser
	 *  \param int: M la taille du vecteur
	 * \param float: K la constante de l'énoncé
	 * \return float*
	 */
		float* initialize(int, float*, int, float) const;

		/*!
	 *  \brief met a jour le vecteur V a chaque étape
	 *
	 *  Methode qui permet de calculer le vecteur Vn en fonction de Vn+1
	 *
	 * \param float**: matrice M1
	 * \param float**: matrice M2
	 * \param float*: vecteur V
	 * \param int: M nombre de lignes de M1 et M2
	 * \param int: N nombre de colonnes de M1 et M2
	 * \param float*: coefficients de la ligne sous la diagonale de M1
	 * \param float*: vecteur k
	 * \param float*: vector1 vecteur temporaire
	 * \param float: K donné
	 * \param float: r donné
	 * \param float: T donné

	 */
		void calculate(float**&  , float**& , float*, int, int, float*, float*, float*, float, float, float);

		//SCHEMA DIFF FINIES

		/*!
	 *  \brief création des matrices pour la methode des différences finies
	 *
	 *  Permet de créer la matrice avec les coefficients correspondants
	 *
	 *  \param float**: matrice nulle
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice carrée
	 * \param int: sigma
	 */
		virtual void computeMatrixDiffFinies(float**&  , float*, float*, int);


		/*!
	 *  \brief création des coefficients de la matrice pour les differences finies
	 *
	 *  Methode qui permet de calculer les coefficients non nuls des matrices de la méthode des différences finies
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice carrée
	 * \param float: mu prix égal a 1

	 */
		void computeCoefDiffFinies(float*, float*, int, float,float);

		/*!
	 *  \brief fonction permettant d'avoir dT correspondant pour la méthode des différences finies
	 * \return float
	 */
		virtual double get_dT_tilde(float sigma) { return dT * 0.5 * sigma * sigma; }

		/*!
	 *  \brief fonction permettant d'initialiser le vecteur V avec la methode diffferences finies
	 * \param int: id de la méthode
	 * \param float*: V le vecteur à initialiser
	 *  \param int: M la taille du vecteur
	 * \param float: K la constante de l'énoncé
	 * \return float*
	 */
		float* initializeDiffFinies(int, float*, int, float) const;

		/*!
	 *  \brief met a jour le vecteur V a chaque étape
	 *
	 *  Methode qui permet de calculer le vecteur Vn en fonction de Vn+1
	 *
	 * \param float**: matrice M1
	 * \param float*: vecteur V
	 * \param int: M nombre de lignes de M1 et M2
	 * \param int: N nombre de colonnes de M1 et M2
	 * \param float*: coefficients de la ligne sous la diagonale de M1
	 * \param float*: vecteur k
	 * \param float*: vector1 vecteur temporaire
	 * \param float: K donné
	 * \param float: r donné
	 * \param float: T donné
	 * \param float: sigma donné

	 */
		void calculateDiffFinies(float**& , float*, int, int, float*, float*, float*, float, float, float,float);






		/*!
	 *  \brief fonction d'affichage du temps et de l'espace
	 *
	 * \param int: M nombre d'interval
	 * \param int: N nombre d'interval
	 * \param float*: V vecteur calculé à chaque étape
	 *
	 *  Affiche le vecteur calculé
	 */
	    void print(int, int, float*) const;

		/*!
	 *  \brief fonction permettant d'avoir dT correspondant
	 * \return float
	 */
		float get_dT() const { return dT; }

		/*!
 *  \brief fonction permettant d'avoir dL correspondant
 * \return float
 */
		float get_dL() const { return dL; }

		/*!
 *  \brief fonction permettant d'avoir l'id de la méthode: 1 pour un call, 0 pour un put
 * \return int
 */
		int get_id() const { return id; }

		

		/*!
	 *  \brief fonction permettant de calculer la puissance d'un  nombre
	 * \param float : a
	 * \param int: b la puissance
	 * \return float
	 */
		float pow(float a, int b) {
			return Resolv::pow(a,b);
		}

		/*!
	 *  \brief fonction permettant d'avoir le maximum entre 2 nombre
	 * \param float : a
	 * \param float: b
	 * \return float
	 */
		float max(float a, float b) const {
			if (a >= b) return a;
			return b;
		}

		/*!
	 *  \brief fonction permettant d'avoir le minimum entre 2 nombre
	 * \param float : a
	 * \param float: b
	 * \return float
	 */
		float min(float a, float b) const {
			if (a <= b) return a;
			return b;
		}

		/*!
	 *  \brief libération de la mémoire prise
	 *
	 *  Methode qui permet de libérer la mémoire utilisée
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param float**: matrice M1
	 * \param float**: matrice M2
	 * \param int: M
	 * \param float*: Vecteur V
	 * \param float*: vecteur k

	 */
		void deleted(float* a, float* b, float* c, float* d, float** M1, float** M2, int M, float* V, float* k) {
			Resolv::deleted(a, b, c, d, M1, M2, M, V, k);

		}

		/*!
	 *  \brief libération de la mémoire prise
	 *
	 *  Methode qui permet de libérer la mémoire utilisée
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param float**: matrice M1
	 * \param float**: matrice M2
	 * \param int: M
	 * \param float*: Vecteur V
	 * \param float*: vecteur k

	 */
	void deletedDiffFinies(float* a, float* b, float** M1, float** M2, int M, float* V, float* k) {
		Resolv::deletedDiffFinies(a, b, M1, M2, M, V, k);
		}

		/*!
	 *  \brief méthode déocmpoisition LU
	 *
	 *  Methode qui permet de mettre en place la décomposition LU pour rédoure Ax=B avec A tridiagonale
	 *
	 * \param float**: matrice A
	 * \param float*: vecteur B
	 * \param int: n

	 */
		float* decompLU(float**& , float*,int);

		/*!
	 *  \brief méthode produit matrice par vecteur
	 *
	 *  Methode qui permet de calculer A*x avec A matrice et x vecteur
	 *
	 * \param float**: matrice A
	 * \param float*: vecteur x
	 * \param int : n nombre de ligne de la matrice
	 * \param int: m nombre de colonne de la matrice

	 */
		float* product(float** & , float*, int, int);

		
	};

	class Call : public ensiie::Resolv {
	private:
		int id;/*!< id de la méthode*/
	public:

		/*!
	 *  \brief Constructeur
	 *
	 *  Constructeur de la classe Call
	 *
	 *  \param float: T le temps
	 *  \param float: L l'espace
	 *  \param int: N nombre d'interval pour l'espace
	 *   \param int M: nombre d'interval pour le temps
	 * \param int: id de la méthode
	 */
		Call(float, float, int, int, int);

		/*!
	 *  \brief Destructeur
	 *
	 *  Destructeur de la classe Call
	 */
		~Call() {};

		/*!
	 *  \brief création des coefficients de la matrice
	 *
	 *  Methode qui permet de calculer les coefficients non nuls des matrices
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice
	 * \param float: sigma donné dans l'énoncé
	 * \param float: r donné dans l'énoncé

	 */
		void computeCoef(float*, float*, float*, float*, int, float, float); 

		/*!
	 *  \brief création des matrices
	 *
	 *  Permet de créer les matrices avec les coefficients correspondants et calculés précédement
	 *
	 *  \param float**: matrice vide
	 * \param float**: matrice vide
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice
	 */
		void computeMatrix(float** , float**, float*, float*, float*, float*, int);

		/*!
	 *  \brief fonction d'affichage du temps et de l'espace
	 *
	 * \param int: M nombre d'interval
	 * \param int: N nombre d'interval
	 * \param float*: V vecteur calculé à chaque étape
	 *
	 *  Affiche le vecteur calculé
	 */
		void print(int, int, float*) const;

		/*!
	 *  \brief fonction permettant d'avoir dT correspondant
	 * \return float
	 */
		float get_dT() const { return dT; }

		/*!
 *  \brief fonction permettant d'avoir dL correspondant
 * \return float
 */
		float get_dL() const { return dL; }

		/*!
 *  \brief fonction permettant d'avoir l'id de la méthode: 1 pour un call, 0 pour un put
 * \return int
 */
		int get_id() const { return id; }

		/*!
	 *  \brief fonction permettant d'initialiser le vecteur V
	 * \param int: id de la méthode
	 * \param float*: V le vecteur à initialiser
	 *  \param int: M la taille du vecteur
	 * \param float: K la constante de l'énoncé
	 * \return float*
	 */
		float* initialize(int, float*, int, float) const;

		/*!
	 *  \brief fonction permettant de calculer la puissance d'un  nombre
	 * \param float : a
	 * \param int: b la puissance
	 * \return float
	 */
		float pow(float a, int b) {
			return Resolv::pow(a, b);
		}

		/*!
	 *  \brief fonction permettant d'avoir le maximum entre 2 nombre
	 * \param float : a
	 * \param float: b
	 * \return float
	 */
		float max(float a, float b) const {
			if (a >= b) return a;
			return b;
		}

		/*!
	 *  \brief fonction permettant d'avoir le minimum entre 2 nombre
	 * \param float : a
	 * \param float: b
	 * \return float
	 */
		float min(float a, float b) const {
			if (a <= b) return a;
			return b;
		}

		/*!
	 *  \brief libération de la mémoire prise
	 *
	 *  Methode qui permet de libérer la mémoire utilisée
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param float**: matrice M1
	 * \param float**: matrice M2
	 * \param int: M
	 * \param float*: Vecteur V
	 * \param float*: vecteur k

	 */
		void deleted(float* a, float* b, float* c, float* d, float** M1, float** M2, int M, float* V, float* k) {
			Resolv::deleted(a, b, c, d, M1, M2, M, V, k);

		}

		/*!
	 *  \brief méthode déocmpoisition LU
	 *
	 *  Methode qui permet de mettre en place la décomposition LU pour rédoure Ax=B avec A tridiagonale
	 *
	 * \param float**: matrice A
	 * \param float*: vecteur B
	 * \param int: n

	 */
		float* decompLU(float**& , float*,int);

		/*!
	 *  \brief méthode produit matrice par vecteur
	 *
	 *  Methode qui permet de calculer A*x avec A matrice et x vecteur
	 *
	 * \param float**: matrice A
	 * \param float*: vecteur x
	 * \param int : n nombre de ligne de la matrice
	 * \param int: m nombre de colonne de la matrice

	 */
		float* product(float**& , float*, int, int);

		/*!
	 *  \brief met a jour le vecteur V a chaque étape
	 *
	 *  Methode qui permet de calculer le vecteur Vn en fonction de Vn+1
	 *
	 * \param float**: matrice M1
	 * \param float**: matrice M2
	 * \param float*: vecteur V
	 * \param int: M nombre de lignes de M1 et M2
	 * \param int: N nombre de colonnes de M1 et M2
	 * \param float*: coefficients de la ligne sous la diagonale de M1
	 * \param float*: vecteur k
	 * \param float*: vector1 vecteur temporaire
	 * \param float: K donné
	 * \param float: r donné
	 * \param float: T donné

	 */
		void calculate(float**& , float**& , float*, int, int, float*, float*, float*, float, float, float);

		//SCHEMA DIFF FINIES

		/*!
	 *  \brief création des matrices pour la methode des différences finies
	 *
	 *  Permet de créer la matrice avec les coefficients correspondants
	 *
	 *  \param float**: matrice nulle
	 *   \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice carrée
	 */
		virtual void computeMatrixDiffFinies(float**, float*, float*, int);


		/*!
	 *  \brief création des coefficients de la matrice pour les differences finies
	 *
	 *  Methode qui permet de calculer les coefficients non nuls des matrices de la méthode des différences finies
	 *
	 *  \param float*: vecteur de coefficients
	 * \param float*: vecteur de coefficients
	 * \param int: M taille de la matrice carrée
	 * \param float: mu 

	 */
		void computeCoefDiffFinies(float*, float*, int, float, float);

		/*!
	 *  \brief fonction permettant d'avoir dT correspondant pour la méthode des différences finies
	 * \return float
	 */
		virtual double  get_dT_tilde(float sigma) { return dT * 0.5 * sigma * sigma; }

		/*!
	 *  \brief fonction permettant d'initialiser le vecteur V avec la methode diffferences finies
	 * \param int: id de la méthode
	 * \param float*: V le vecteur à initialiser
	 *  \param int: M la taille du vecteur
	 * \param float: K la constante de l'énoncé
	 * \return float*
	 */
		float* initializeDiffFinies(int, float*, int, float) const;

		/*!
	 *  \brief met a jour le vecteur V a chaque étape
	 *
	 *  Methode qui permet de calculer le vecteur Vn en fonction de Vn+1
	 *
	 * \param float**: matrice M1
	 * \param float*: vecteur V
	 * \param int: M nombre de lignes de M1 et M2
	 * \param int: N nombre de colonnes de M1 et M2
	 * \param float*: coefficients de la ligne sous la diagonale de M1
	 * \param float*: vecteur k
	 * \param float*: vector1 vecteur temporaire
	 * \param float: K donné
	 * \param float: r donné
	 * \param float: T donné
	 * \param float: sigma donné

	 */
		void calculateDiffFinies(float**& , float*, int, int, float*, float*, float*, float, float, float, float);


	};
}


#endif