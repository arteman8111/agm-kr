const render = (arr, label, index) => {
    const app = document.querySelector('#app')
    const table = document.createElement('table')
    table.className = "data"
    const header = document.createElement('tr')
    const td = (element) => `<td class="data__element">${element}</td>`
    header.insertAdjacentHTML('beforeend',td(''))
    for(let i = index; i < arr[0].length + index; i++){
        header.insertAdjacentHTML('beforeend', td(i))
    }
    for (let i = 0; i < arr.length; i++) {
        const tr = document.createElement('tr');
        tr.insertAdjacentHTML('beforeend',td(label[i]))
        arr[i].forEach(element => tr.insertAdjacentHTML('beforeend', td(math.format(element,3))))
        table.insertAdjacentElement('beforeend', tr)
    }
    table.insertAdjacentElement('afterbegin', header)
    app.insertAdjacentElement('beforeend', table)
}

const id = () => {
        
}

export { render, id }