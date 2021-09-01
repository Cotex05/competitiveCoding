#include <iostream>
#include <math.h>
using namespace std;

template <class T>
struct Node {
    T value;
    Node *left;
    Node *right;

    Node(T val) {
        this->value = val;
    }

    Node(T val, Node<T> left, Node<T> right) {
        this->value = val;
        this->left = left;
        this->right = right;
    }
};

template <class T>
class BST {

    private:
    Node<T> *root;

    void addHelper(Node<T> *root, T val) {
        if (root->value > val) {
            if (!root->left) {
                root->left = new Node<T>(val);
            } else {
                addHelper(root->left, val);
            }
        } else {
            if (!root->right) {
                root->right = new Node<T>(val);
            } else {
                addHelper(root->right, val);
            }
        }
    }

    void printInorderHelper(Node<T> *root) {
        if (!root) return;
        printInorderHelper(root->left);
        cout<<root->value<<' ';
        printInorderHelper(root->right);
    }

    void printPreorderHelper(Node<T> *root) {
        if (!root) return;
        cout<<root->value<<' ';
        printPreorderHelper(root->left);
        printPreorderHelper(root->right);
    }
    
    void printPostorderHelper(Node<T> *root) {
        if (!root) return;
        printPostorderHelper(root->left);
        printPostorderHelper(root->right);
        cout<<root->value<<' ';
    }
    
    void printLevelorderHelper(Node<T> *root)
    {
        T h = heightHelper(root);
        T i;
        for (i = 1; i <= h; i++)
            printCurrentLevel(root, i);
    }
     
    /* Print nodes at a current level */
    void printCurrentLevel(Node<T> *root, T level)
    {
        if (root == NULL)
            return;
        if (level == 1)
            cout << root->value << " ";
        else if (level > 1)
        {
            printCurrentLevel(root->left, level-1);
            printCurrentLevel(root->right, level-1);
        }
    }

    int nodesCountHelper(Node<T> *root) {
        if (!root) return 0;
        else return 1 + nodesCountHelper(root->left) + nodesCountHelper(root->right);
    }

    int heightHelper(Node<T> *root) {
        if (!root) return 0;
        else return 1 + max(heightHelper(root->left), heightHelper(root->right));
    }

    void printMaxPathHelper(Node<T> *root) {
        if (!root) return;
        cout<<root->value<<' ';
        if (heightHelper(root->left) > heightHelper(root->right)) {
            printMaxPathHelper(root->left);
        } else {
            printMaxPathHelper(root->right);
        }
    }

    bool deleteValueHelper(Node<T>* parent, Node<T>* current, T value) {
        if (!current) return false;
        if (current->value == value) {
            if (current->left == NULL || current->right == NULL) {
                Node<T>* temp = current->left;
                if (current->right) temp = current->right;
                if (parent) {
                    if (parent->left == current) {
                        parent->left = temp;
                    } else {
                        parent->right = temp;
                    }
                } else {
                    this->root = temp;
                }
            } else {
                Node<T>* validSubs = current->right;
                while (validSubs->left) {
                    validSubs = validSubs->left;
                }
                T temp = current->value;
                current->value = validSubs->value;
                validSubs->value = temp;
                return deleteValueHelper(current, current->right, temp);
            }
            delete current;
            return true;
        }
        return deleteValueHelper(current, current->left, value) ||
               deleteValueHelper(current, current->right, value);
    }

    public:
    void add(T val) {
        if (root) {
            this->addHelper(root, val);
        } else {
            root = new Node<T>(val);
        }
    }

    void printInorder() {
        printInorderHelper(this->root); 
    }
    
    void printPreorder() {
        printPreorderHelper(this->root);
    }
    
    void printPostorder() {
        printPostorderHelper(this->root);
    }
    
    void printLevelorder(){
        printLevelorderHelper(this->root);
    }

    int nodesCount() {
        return nodesCountHelper(root);
    }

    int height() {
        return heightHelper(this->root);
    }

    void printMaxPath() {
        printMaxPathHelper(this->root);
    }

    bool deleteValue(T value) {
        return this->deleteValueHelper(NULL, this->root, value);
    }
};

int main() {
    BST<int> *bst = new BST<int>();
    bst->add(10);
    bst->add(15);
    bst->add(8);
    bst->add(5);
    bst->add(12);
    bst->add(33);
    bst->add(7);
    bst->printInorder();
    cout<<endl;
    bst->printPreorder();
    cout<<endl;
    bst->printPostorder();
    cout<<endl;
    bst->printLevelorder();
    cout<<endl;
    return 0;
}
