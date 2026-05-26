Foam::string commonSubstring(const Foam::string& a, const Foam::string& b)
{
    Foam::string common;

    for (size_t i = 0; i < a.size(); ++i)
    {
        for (size_t j = i; j < a.size(); ++j)
        {
            Foam::string sub = a.substr(i, j - i + 1);
            if 
            (
                b.find(sub) 
             != Foam::string::npos && sub.size() > common.size()
            )
            {
                common = sub;
            }
        }
    }

    return common;
}
