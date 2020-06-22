#include "Matrix.hpp"

// Constructor for Any Matrix
Matrix::Matrix(unsigned rowSize, unsigned colSize, double initial) {
    rows = rowSize;
    cols = colSize;
    data.resize(rowSize);
    for (unsigned i = 0; i < data.size(); i++) {
        data[i].resize(colSize, initial);
    }
}

// Constructor for Any Matrix
Matrix::Matrix(unsigned rowSize, unsigned colSize) : Matrix(rowSize, colSize, 0.0f) {}

// Constructor for Given Matrix
Matrix::Matrix(const char *fileName) {
    std::ifstream file_A(fileName);
    int colSize = 0;
    int rowSize = 0;
    std::string line_A;
    int idx = 0;
    double element_A;
    double *vector_A = nullptr;

    if (fileName == "Error") {
        std::cout << "The dimentions of Matrix are wrong for this operation!\n";
    } else {
        if (file_A.is_open() && file_A.good()) {
            // cout << "File A.txt is open. \n";
            while (getline(file_A, line_A)) {
                rowSize += 1;
                std::stringstream stream_A(line_A);
                colSize = 0;
                while (1) {
                    stream_A >> element_A;
                    if (!stream_A)
                        break;
                    colSize += 1;
                    double *tempArr = new double[idx + 1];
                    std::copy(vector_A, vector_A + idx, tempArr);
                    tempArr[idx] = element_A;
                    vector_A = tempArr;
                    idx += 1;
                }
            }
        } else {
            std::cout << " Failed to open.\n";
        }
    }
    int j;
    idx = 0;
    data.resize(rowSize);
    for (unsigned i = 0; i < data.size(); i++) {
        data[i].resize(colSize);
    }
    for (int i = 0; i < rowSize; i++) {
        for (j = 0; j < colSize; j++) {
            this->data[i][j] = vector_A[idx];
            idx++;
        }
    }
    cols = colSize;
    rows = rowSize;
    delete[] vector_A;
}

// Copy Constructor
Matrix::Matrix(const Matrix &B) {
    this->cols = B.getCols();
    this->rows = B.getRows();
    this->data = B.data;
}

// From Array Constructor
Matrix::Matrix(vector<double> array) {
    this->cols = 1;
    this->rows = array.size();
    this->data.resize(rows);
    for (unsigned i = 0; i < data.size(); i++) {
        data[i].resize(cols, array[i]);
    }
}

Matrix::Matrix() {
    this->rows = 0;
    this->cols = 0;
}

// Addition of Two Matrices
Matrix Matrix::operator+(Matrix &B) {
    Matrix sum(rows, cols, 0.0);
    //sum.print();
    if ((this->rows != B.getRows()) || (this->cols != B.getCols())) {
        return "Error";
    }
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            sum(i, j) = data[i][j] + B(i, j);
        }
    }
    return sum;
}

// Subtraction of Two Matrices
Matrix Matrix::operator-(Matrix &B) {
    Matrix diff(rows, cols, 0.0);
    if ((this->rows != B.getRows()) || (this->cols != B.getCols())) {
        return "Error";
    }
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            diff(i, j) = this->data[i][j] - B(i, j);
        }
    }

    return diff;
}

// Multiplication of Two Matrices
Matrix Matrix::operator*(Matrix &B) {
    Matrix multip(this->rows, B.getCols(), 0.0);
    if (this->cols == B.getRows()) {
        unsigned i, j, k;
        double temp = 0.0;
        for (i = 0; i < this->rows; i++) {
            for (j = 0; j < B.getCols(); j++) {
                temp = 0.0;
                for (k = 0; k < this->cols; k++) {
                    temp += data[i][k] * B(k, j);
                }
                multip(i, j) = temp;
            }
        }
        return multip;
    } else {
        return "Error";
    }
}

// Addition of Two Matrices
Matrix Matrix::operator+=(Matrix &B) {
    //sum.print();
    if ((this->rows != B.getRows()) || (this->cols != B.getCols())) {
        return "Error";
    }
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            this->data[i][j] = this->data[i][j] + B(i, j);
        }
    }
    return *this;
}

// Subtraction of Two Matrices
Matrix Matrix::operator-=(Matrix &B) {
    if ((this->rows != B.getRows()) || (this->cols != B.getCols())) {
        return "Error";
    }
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            this->data[i][j] = this->data[i][j] - B(i, j);
        }
    }

    return *this;
}

// Multiplication of Two Matrices
Matrix Matrix::operator*=(Matrix &B) {
    Matrix multip(this->rows, B.getCols(), 0.0);
    if (this->cols == B.getRows()) {
        unsigned i, j, k;
        double temp = 0.0;
        for (i = 0; i < this->rows; i++) {
            for (j = 0; j < B.getCols(); j++) {
                temp = 0.0;
                for (k = 0; k < this->cols; k++) {
                    temp += data[i][k] * B(k, j);
                }
                multip(i, j) = temp;
            }
        }
        *this = multip;
        return *this;
    } else {
        return "Error";
    }
}

// Scalar Addition
Matrix Matrix::operator+(double scalar) {
    Matrix result(rows, cols, 0.0);
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            result(i, j) = this->data[i][j] + scalar;
        }
    }
    return result;
}

// Scalar Subraction
Matrix Matrix::operator-(double scalar) {
    Matrix result(rows, cols, 0.0);
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            result(i, j) = this->data[i][j] - scalar;
        }
    }
    return result;
}

// Scalar Multiplication
Matrix Matrix::operator*(double scalar) {
    Matrix result(rows, cols, 0.0);
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            result(i, j) = this->data[i][j] * scalar;
        }
    }
    return result;
}

// Scalar Division
Matrix Matrix::operator/(double scalar) {
    Matrix result(rows, cols, 0.0);
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            result(i, j) = this->data[i][j] / scalar;
        }
    }
    return result;
}

// Scalar Addition
Matrix Matrix::operator+=(double scalar) {
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            this->data[i][j] = this->data[i][j] + scalar;
        }
    }
    return *this;
}

// Scalar Subraction
Matrix Matrix::operator-=(double scalar) {
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            this->data[i][j] = this->data[i][j] - scalar;
        }
    }
    return *this;
}

// Scalar Multiplication
Matrix Matrix::operator*=(double scalar) {
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            this->data[i][j] = this->data[i][j] * scalar;
        }
    }
    return *this;
}

// Scalar Division
Matrix Matrix::operator/=(double scalar) {
    unsigned i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            this->data[i][j] = this->data[i][j] / scalar;
        }
    }
    return *this;
}

// Returns value of given location when asked in the form A(x,y)
double &Matrix::operator()(const unsigned &rowNo, const unsigned &colNo) {
    return this->data[rowNo][colNo];
}

// No brainer - returns row #
unsigned Matrix::getRows() const {
    return this->rows;
}

// returns col #
unsigned Matrix::getCols() const {
    return this->cols;
}

vector<double> Matrix::toArray() {
    vector<double> array;
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            array.push_back(data[i][j]);
        }
    }
    return array;
}

// Take any given matrices transpose and returns another matrix
Matrix Matrix::transpose() {
    Matrix Transpose(rows, cols, 0.0);
    for (unsigned i = 0; i < cols; i++) {
        for (unsigned j = 0; j < rows; j++) {
            Transpose(i, j) = this->data[j][i];
        }
    }
    return Transpose;
}

// Prints the matrix beautifully
void Matrix::print() const {
    if (this->rows == 0) std::cout << "[Empty Matrix]\n";
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            std::cout << "[" << data[i][j] << "] ";
        }
        std::cout << '\n';
    }
}
// Prints the matrix beautifully
std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
    if (matrix.rows == 0) os << "[Empty Matrix]\n";
    for (unsigned i = 0; i < matrix.rows; i++) {
        for (unsigned j = 0; j < matrix.cols; j++) {
            os << "[" << matrix.data[i][j] << "] ";
        }
        os << '\n';
    }
    return os;
}

// Returns 3 values
//First: Eigen Vector
//Second: Eigen Value
//Third: Flag
tuple<Matrix, double, int> Matrix::powerIter(unsigned rowNum, double tolerance) {
    // Picks a classic X vector
    Matrix X(rowNum, 1, 1.0);
    // Initiates X vector with values 1,2,3,4
    for (unsigned i = 1; i <= rowNum; i++) {
        X(i - 1, 0) = i;
    }
    int errorCode = 0;
    double difference = 1.0;  // Initiall value greater than tolerance
    unsigned j = 0;
    unsigned location;
    // Defined to find the value between last two eigen values
    vector<double> eigen;
    double eigenvalue = 0.0;
    eigen.push_back(0.0);

    while (abs(difference) > tolerance) {  // breaks out when reached tolerance
        j++;
        // Normalize X vector with infinite norm
        for (int i = 0; i < rowNum; ++i) {
            eigenvalue = X(0, 0);
            if (abs(X(i, 0)) >= abs(eigenvalue)) {
                // Take the value of the infinite norm as your eigenvalue
                eigenvalue = X(i, 0);
                location = i;
            }
        }
        if (j >= 5e5) {
            std::cout << "Oops, that was a nasty complex number wasn't it?\n";
            std::cout << "ERROR! Returning code black, code black!";
            errorCode = -1;
            return std::make_tuple(X, 0.0, errorCode);
        }
        eigen.push_back(eigenvalue);
        difference = eigen[j] - eigen[j - 1];
        // Normalize X vector with its infinite norm
        X = X / eigenvalue;

        // Multiply The matrix with X vector
        X = (*this) * X;
    }

    // Take the X vector and what you've found is an eigenvector!
    X = X / eigenvalue;
    return std::make_tuple(X, eigenvalue, errorCode);
}

Matrix Matrix::deflation(Matrix &X, double &eigenvalue) {
    // Deflation formula exactly applied
    double denominator = eigenvalue / (X.transpose() * X)(0, 0);
    Matrix Xtrans = X.transpose();
    Matrix RHS = (X * Xtrans);
    Matrix RHS2 = RHS * denominator;
    Matrix A2 = *this - RHS2;
    return A2;
}

// Randomize all elements of matrix (range [-1,1])
void Matrix::randomize() {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> dist(-1.0f, 1.0f);
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            double value = this->data[i][j];
            this->data[i][j] = dist(rng);
        }
    }
}

//Activation Function
Matrix Matrix::sigmod() {
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            double value = this->data[i][j];
            this->data[i][j] = 1 / (1 + std::exp(-value));
        }
    }

    return *this;
}

// Gradient function
Matrix Matrix::dSigmod() {
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            double value = this->data[i][j];
            this->data[i][j] = 1 - (value * value);
        }
    }

    return *this;
}

// Apply a function in matrix
Matrix Matrix::map(functioncall func) {
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            double value = this->data[i][j];
            this->data[i][j] = func(value);
        }
    }

    return *this;
}