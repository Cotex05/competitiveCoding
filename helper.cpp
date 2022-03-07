#include <bits/stdc++.h>
#define ll long long int
#define tc    \
    ll t;     \
    cin >> t; \
    while (t--)
#define all(x) x.begin(), x.end()
#define deb(x) cout << #x << " " << x << "\n";
const int INF = INT_MAX;
#define INFL LLONG_MAX
#define pb push_back
#define eb emplace_back
#define sum(a) (accumulate((a).begin(), (a).end(), 0ll))
#define F(i, a, n) for (int i = a; i < n; i++)
#define mod 1000000007
using namespace std;

// Useful Headers
// #include <boost/multiprecision/cpp_int.hpp>
// using boost::multiprecision::cpp_int;

template <typename... T>
void read(T &...args)
{
    ((cin >> args), ...);
}

bool cmp(pair<int, int> a, pair<int, int> b)
{
    return a.second < b.second; // ascending order by second element
}

bool isPalindrome(string s)
{
    int i = 0, j = s.size() - 1;
    while (i < j)
    {
        if (s[i] != s[j])
        {
            return false;
        }
        i++;
        j--;
    }
    return true;
}

bool isPowerOfTwo(ll x)
{
    return x && (!(x & (x - 1)));
}

bool isPrime(ll x)
{
    for (ll i = 2; i <= x / 2; i++)
    {
        if (x % i == 0)
        {
            return false;
        }
    }
    return true;
}

// count all primes below n
int countPrimes(int n)
{
    if (n < 1)
    {
        return 0;
    }
    vector<bool> prime(n, true);
    prime[0] = false, prime[1] = false;
    for (int i = 0; i < sqrt(n); ++i)
    {
        if (prime[i])
        {
            for (int j = i * i; j < n; j += i)
            {
                prime[j] = false;
            }
        }
    }
    return count(prime.begin(), prime.end(), true);
}

vector<ll> factors(ll x)
{
    vector<ll> result;
    ll i = 1;
    while (i * i <= x)
    {
        if (x % i == 0)
        {
            result.push_back(i);
            if (x / i != i)
            {
                result.push_back(x / i);
            }
        }
        i++;
    }
    return result;
}

vector<ll> findPrimeFactors(ll n)
{
    vector<ll> result;
    for (ll i = 2; i <= n; i++)
    {
        while (n % i == 0)
        {
            result.push_back(i);
            n /= i;
        }
    }
    return result;
}

ll countPrimeFactors(ll N)
{
    ll count = 0;
    for (ll i = 2; i * i <= N; ++i)
    {
        if (N % i == 0)
        {
            count++;
            while (N % i == 0)
            {
                N = N / i;
            }
        }
    }
    if (N > 1)
        count++;
    return count;
}

ll power(ll x, ll y)
{
    ll result = 1;
    while (y > 0)
    {
        if (y % 2 == 0) // y is even
        {
            x = x * x;
            y = y / 2;
        }
        else // y isn't even
        {
            result = result * x;
            y = y - 1;
        }
    }
    return result;
}

ll factorial(ll n)
{
    if (n < 0)
        return (-1); /*Wrong value*/
    if (n == 0)
        return (1); /*Terminating condition*/
    else
    {
        return (n * factorial(n - 1));
    }
}

ll lcm(ll x, ll y)
{
    return x * y / __gcd(x, y);
}

long long subarrayXor(int arr[], int n)
{
    long long ans = 0;

    int *xorArr = new int[n];

    unordered_map<int, int> mp;

    xorArr[0] = arr[0];

    for (int i = 1; i < n; i++)
        xorArr[i] = xorArr[i - 1] ^ arr[i];

    for (int i = 0; i < n; i++)
    {

        int tmp = 0 ^ xorArr[i];

        ans = ans + ((long long)mp[tmp]);

        if (xorArr[i] == 0)
            ans++;

        mp[xorArr[i]]++;
    }

    return ans;
}

int max_xor(int arr[], int n)
{
    int maxx = 0, mask = 0;

    set<int> se;

    for (int i = 30; i >= 0; i--)
    {

        mask |= (1 << i);

        for (int i = 0; i < n; ++i)
        {

            se.insert(arr[i] & mask);
        }

        int newMaxx = maxx | (1 << i);

        for (int prefix : se)
        {

            if (se.count(newMaxx ^ prefix))
            {
                maxx = newMaxx;
                break;
            }
        }

        se.clear();
    }

    return maxx;
}

int hamming(int x, int y)
{
    return __builtin_popcount(x ^ y);
}

// A simple method to evaluate Euler Totient Function
// Properties
// -> phi(p)=p-1 , p is prime.
// -> phi(a*b) = phi(a)*phi(b) , a and b are primes.
// -> IF gcd(a, b) = 1, then phi(a*b) = phi(a)*phi(b)*gcd(a,b)/phi(gcd(a,b)),  a and b are primes.
// -> a^(p-1) ≡ 1 (mod p), fermat's little theorem.
// -> a^Φ(n) ≡ 1 (mod n)
//  O(nlogn)
ll phi(unsigned ll n)
{
    unsigned ll result = 1;
    for (ll i = 2; i < n; i++)
        if (__gcd(i, n) == 1)
            result++;
    return result;
}

// O(n)
ll phi(ll n)
{
    ll result = n;
    for (ll p = 2; p * p <= n; ++p)
    {
        if (n % p == 0)
        {
            while (n % p == 0)
                n /= p;

            result -= result / p;
        }
    }
    if (n > 1)
        result -= result / n;

    return result;
}

// function generate all prime number less then N in O(n)
const long long MAX_SIZE = 1000001;

vector<ll> isprime(MAX_SIZE, true);
vector<ll> prime;
vector<ll> SPF(MAX_SIZE);

void manipulatedSeive(int N)
{
    isprime[0] = isprime[1] = false;

    for (long long int i = 2; i < N; i++)
    {
        if (isprime[i])
        {
            prime.push_back(i);
            SPF[i] = i;
        }
        for (long long int j = 0;
             j < (int)prime.size() &&
             i * prime[j] < N && prime[j] <= SPF[i];
             j++)
        {
            isprime[i * prime[j]] = false;

            SPF[i * prime[j]] = prime[j];
        }
    }
}

// Geometry
struct Point
{
    int x;
    int y;
};

// Given three collinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(Point p, Point q, Point r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;
    return false;
}

int orientation(Point p, Point q, Point r)
{
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0)
        return 0;             // collinear
    return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// The function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(Point p1, Point q1, Point p2, Point q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1))
        return true;

    // p1, q1 and p2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1))
        return true;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2))
        return true;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2))
        return true;

    return false; // Doesn't fall in any of the above cases
}

// Graph

// Taking input i.e Adjacency List representaion
void adjInGraph(ll n, ll m)
{
    cin >> n >> m;
    vector<ll> adj[n + 1];
    for (ll i = 0; i < m; i++)
    {
        ll u, v;
        cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
}

// BFS (start node, adj. List, total nodes/vertices)
void BFS(ll s, vector<ll> adj[], ll n)
{
    bool *visited = new bool[n];
    for (ll i = 0; i < n; i++)
        visited[i] = false;

    list<ll> queue;
    visited[s] = true;
    queue.push_back(s);

    while (!queue.empty())
    {
        s = queue.front();
        cout << s << " ";
        queue.pop_front();

        for (auto i : adj[s])
        {
            if (!visited[i])
            {
                visited[i] = true;
                queue.push_back(i);
            }
        }
    }
}

// DFS
void DFS(ll s, vector<ll> adj[], ll n)
{
    visited[s] = true;
    cout << s << " ";
    for (auto it : adj[s])
        if (!visited[it])
            DFS(it, adj, n);
}

void dfs(ll s)
{
    visited[s] = true;
    for (ll child : adj[s])
    {
        if (visited[child])
            continue;
        dfs(child);
    }
}

// Number of components in graph
void components(ll n)
{
    ll count = 0;
    for (ll i = 0; i < n; i++)
    {
        if (visited[i])
            continue;
        dfs(i);
        count++;
    }
}

// Detect Cycle in a undirected graph using DFS
class Solution
{
public:
    bool checkForCycle(ll node, ll parent, vector<bool> &vis, vector<ll> adj[])
    {
        vis[node] = true;
        for (auto it : adj[node])
        {
            if (!vis[it])
            {
                if (checkForCycle(it, node, vis, adj))
                    return true;
            }
            else if (it != parent)
                return true;
        }
        return false;
    }
    bool isCycle(ll n, vector<ll> adj[])
    {
        vector<bool> vis(n + 1, false);
        for (ll i = 1; i <= n; i++)
        {
            if (!vis[i])
            {
                if (checkForCycle(i, -1, vis, adj))
                    return true;
            }
        }
        return false;
    }
}

// Detect Cycle using BFS
class Solution
{
public:
    bool checkForCycle(ll s, ll n, vector<bool> &vis, vector<ll> adj[])
    {
        // queue for BFS...
        queue<pair<ll, ll>> q;
        vis[s] = true;
        q.push({s, -1});

        while (!q.empty())
        {
            ll node = q.front().first;
            ll par = q.front().second;
            q.pop();

            for (auto it : adj[node])
            {
                if (!vis[it])
                {
                    vis[it] = true;
                    q.push({it, node});
                }
                else if (node != par)
                {
                    return true;
                }
            }
        }
        return false;
    }
    bool isCycle(ll n, vector<ll> adj[])
    {
        vector<bool> vis(n + 1, false);
        for (ll i = 1; i <= n; i++)
        {
            if (!vis[i])
            {
                if (checkForCycle(i, n, vis, adj))
                {
                    return true;
                }
            }
        }
        return false;
    }
}

class findLCA
{
    int N = 1000;

public:
    vector<int> gr[N];
    int dep[N], Par[N];

    void dfs(int cur, int par)
    {
        Par[cur] = par;
        dep[cur] = dep[par] + 1;
        for (auto x : gr[cur])
        {
            if (x != par)
            {
                dfs(x, cur);
            }
        }
    }

    int LCA(int u, int v)
    {
        if (u == v)
            return u;

        if (dep[u] < dep[v])
            swap(u, v);
        // depth of u is more than depth of v

        int diff = dep[u] - dep[v];

        // depth of both nodes same
        while (diff--)
        {
            u = Par[u];
        }

        // until they are equal nodes keep climbing
        while (u != v)
        {
            u = Par[u];
            v = Par[v];
        }

        return u;
    }

    // dfs(1, 0);

    // LCA(9, 6)
}

// Check Bipartite graph using DFS
bool
bipartiteDfs(int node, vector<int> adj[], int color[])
{
    for (auto it : adj[node])
    {
        if (color[it] == -1)
        {
            color[it] = 1 - color[node];
            if (!bipartiteDfs(it, adj, color))
            {
                return false;
            }
        }
        else if (color[it] == color[node])
            return false;
    }
    return true;
}

bool checkBipartite(vector<int> adj[], int n)
{
    int color[n];
    memset(color, -1, sizeof color);
    for (int i = 0; i < n; i++)
    {
        if (color[i] == -1)
        {
            color[i] = 1;
            if (!bipartiteDfs(i, adj, color))
            {
                return false;
            }
        }
    }
    return true;
}

// DAG -> Directed acyclic graph
// Topological Sort/Order using DFS
void findTopoSort(int node, vector<int> &vis, stack<int> &st, vector<int> adj[])
{
    vis[node] = 1;
    for (auto it : adj[node])
    {
        if (!vis[it])
        {
            findTopoSort(it, vis, st, adj);
        }
    }
    st.push(node);
}

vector<int> topoSort(int V, vector<int> adj[])
{
    // code here
    stack<int> st;
    vector<int> vis(V, 0);
    for (int i = 0; i < V; i++)
    {
        if (!vis[i])
        {
            findTopoSort(i, vis, st, adj);
        }
    }
    vector<int> ans;
    while (!st.empty())
    {
        int top = st.top();
        ans.push_back(top);
        st.pop();
    }
    return ans;
}

// Topological Sort using BFS
vector<int> topologicalSort(int n, vector<int> adj[])
{
    queue<int> q;
    vector<int> inDegree(n, 0);
    for (int i = 0; i < n; i++)
    {
        for (auto it : adj[i])
        {
            inDegree[it]++;
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (inDegree[i] == 0)
        {
            q.push(i);
        }
    }
    vector<int> topo;
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        topo.push_back(node);
        for (auto it : adj[node])
        {
            inDegree[it]--;
            if (inDegree[it] == 0)
            {
                q.push(it);
            }
        }
    }
    return topo;
}

// Cycle detection in directed graph using topo-sort
// if counts of topo array is equal to total nodes then no cycle is present.
//  i.e. if(topo.size() == n) return false;

// Shortest Path in Undirected Graph with Unit Weights using BFS

vector<int> shortestDist(int n, vector<int> adj[], int src)
{
    vector<int> dist(n, INT_MAX);
    queue<int> q;

    dist[src] = 0;
    q.push(src);

    while (!q.empty())
    {
        int node = q.front();
        q.pop();

        for (auto it : adj[node])
        {
            if (dist[node] + 1 < dist[it])
            {
                dist[it] = dist[node] + 1;
                q.push(it);
            }
        }
    }
    return dist;
}

// DSU
class Graph
{

    list<pair<int, int>> edge_list;
    int V;

public:
    Graph(int V)
    {
        this->V = V;
    }

    void addEdge(int u, int v)
    {
        edge_list.push_back(make_pair(u, v));
    }

    int findSet(int i, int parent[])
    {
        if (parent[i] == -1)
        {
            return i;
        }
        int p = findSet(parent[i], parent);
        // compress the path for next time
        parent[i] = p;
        return p;
    }
    void union_set(int u, int v, int rank[], int parent[])
    {

        int p1 = findSet(u, parent);
        int p2 = findSet(v, parent);

        if (p1 != p2)
        {
            if (rank[p1] < rank[p2])
            {
                parent[p1] = p2;
                rank[p2] += rank[p1];
            }
            else
            {
                parent[p2] = p1;
                rank[p1] += rank[p2];
            }
        }
        // print parent array after every union
        //  cout<<"Parent";
        // for (int i = 0; i < V; i++)
        // {
        //     cout<<parent[i]<<" ";
        // }
        // cout<<endl;
        // cout<<"Rank";
        // for (int j = 0; j < V; j++)
        // {
        //     cout<<rank[j]<<" ";
        // }
        // cout<<endl;
    }

    bool contains_cycle()
    {
        int *parent = new int[V];
        int *rank = new int[V];
        for (int i = 0; i < V; i++)
        {
            parent[i] = -1;
            rank[i] = 1;
        }

        for (auto edge : edge_list)
        {
            int i = edge.first;
            int j = edge.second;

            int p1 = findSet(i, parent);
            int p2 = findSet(j, parent);
            // cout<<i<<"-"<<p1<<" and "<<j<<"-"<<p2<<endl;

            if (p1 != p2)
            {
                union_set(p1, p2, rank, parent);
            }
            else
            {
                // belong to same set
                //  cout<<"Parent";
                //  for(int i=0;i<V;i++){
                //  	cout<<parent[i]<<" ";
                //  }
                return true;
            }
        }
        return false;
        delete[] parent;
    }
};

// Dijkstra Algorithm
// shortest distance of all the vertex's from the source vertex S
vector<int> dijkstra(int V, vector<vector<int>> adj[], int S)
{
    // Code here
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    vector<int> distTo(V, INT_MAX);

    distTo[S] = 0;
    pq.push(make_pair(0, S)); // weight, node

    while (!pq.empty())
    {
        int dist = pq.top().first;
        int prev = pq.top().second;
        pq.pop();

        for (auto ele : adj[prev])
        {
            int next = ele[0];
            int nextDist = ele[1];

            if (distTo[next] > distTo[prev] + nextDist)
            {
                distTo[next] = distTo[prev] + nextDist;
                pq.push(make_pair(distTo[next], next));
            }
        }
    }
    return distTo;
}

// Bellman Ford's algo with negative edges.
vector<int> bellman_ford(int n, vector<vector<int>> adj, int src)
{
    // Code here
    vector<int> dist(n, 1e8);

    dist[src] = 0;

    for (int i = 0; i < n; i++)
    {
        for (auto &it : adj)
        {
            if (dist[it[0]] + it[2] < dist[it[1]])
            {
                dist[it[1]] = dist[it[0]] + it[2];
            }
        }
    }
    return dist;
}
// Prim's Algorithm for MST
vector<int> prims(int n, vector<pair<int, int>> adj[])
{
    vector<int> parent(n, -1);
    int key[n] = {INT_MAX};
    bool mstSet[n] = {false};

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

    key[0] = 0;
    parent[0] = -1;
    pq.push(make_pair(0, 0));

    while (!pq.empty())
    {
        int u = pq.top().second;
        pq.pop();

        mstSet[u] = true;

        for (auto it : adj[u])
        {
            int v = it.first;
            int weight = it.second;
            if (mstSet[v] == false && weight > key[v])
            {
                parent[v] = u;
                key[v] = weight;
                pq.push(make_pair(key[v], v));
            }
        }
    }
    return parent;
}

// Floyd Warshall
vector<vector<int>> floyd_warshall(vector<vector<int>> adj, int n)
{
    vector<vector<int>> dist(adj);

    for (int k = 0; k < n; k++)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (dist[i][j] > dist[i][k] + dist[k][j])
                {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
    return dist;
}

// Travelling salesman problem
int tsp(vector<vector<int>> &cost, int i, int depth, vector<int> &vis)
{
    int n = cost.size();
    if (depth == n)
    {
        return cost[i][0];
    }
    int ans = INT_MAX;
    vis[i] = 1;
    for (int j = 0; j < n; j++)
    {
        if (!vis[j])
        {
            vis[j] = 1;
            ans = min(ans, cost[i][j] + tsp(cost, j, depth + 1, vis));
            vis[j] = 0;
        }
    }
    vis[i] = 0;
    return ans;
}

int total_cost(vector<vector<int>> cost)
{
    vector<int> vis(cost.size(), 0);
    return tsp(cost, 0, 1, vis);
}