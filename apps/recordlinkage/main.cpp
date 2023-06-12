int main(int argc, char* argv[]) {
    for (int i = 0; i < 10; i++) {
        auto array = new bool[900]();
        for (int ii=0; ii<900; ii++) {
            if (array[ii]) {
                cout << "not zero initialised" << endl;
                return 1;
            }
        }
        delete[] array;
    }
}
