#ifndef _my_parameter_handler_h_
#define _my_parameter_handler_h_

namespace MyParameterHandler
{

    template<typename ItemType>
    std::vector<ItemType> get_vector(ParameterHandler &prm, std::string parameter_name)
    {
        std::vector<std::string> strings = Utilities::split_string_list(prm.get(parameter_name));
        std::vector<ItemType> items;
        for (auto &string : strings) 
        {
            std::stringstream parser(string);
            ItemType item;
            parser >> item;
            items.push_back(item);
        }
        return items;
    } 

}

#endif
